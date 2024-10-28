
def openbis_login(url, username, s3_config_path=None):
    from getpass import getpass

    if s3_config_path:
        from ob.OpenbisAixTended import OpenbisWithS3
        o = OpenbisWithS3(url=url, verify_certificates=True, 
                          s3_config_path=s3_config_path, standardize_filenames=True)
    else:
        from pybis import Openbis
        o = Openbis(url)

    if not o.is_session_active():
        o.login(username, getpass('Enter openBIS password: '), save_token=True) # save the session token in ~/.pybis/example.com.token
    return o

def format_json_string(json_string):
    json_string = json_string.replace('\n', '<br>')
    result = []
    for index, char in enumerate(json_string):
        if char == " " and (index == 0 or json_string[index - 1] != ":"):
            result.append("&nbsp;&nbsp;")
        else:
            result.append(char)

    json_string = "".join(result)
    return json_string

def flatten_lammps_cdict(cdict):
        flat = {}
        for k, v in cdict.items():
            if k != '@context':
                if isinstance(v, dict):
                    if 'label' in v.keys():
                        flat[k] = v['label']
                    else:
                        flat = flat | flatten_lammps_cdict(v)
                elif k == 'inputs' or k == 'outputs' or k == 'job_details':
                    for i in v:
                        flat[i['label']] = i['value']
                elif k == 'software':
                    flat[k] = v[0]['label']
                else:
                    flat[k] = v
        return flat

def flatten_crystal_cdict(cdict):
        flat = {}
        for k, v in cdict.items():
            if k != '@context':
                if isinstance(v, dict):
                    if 'label' in v.keys():
                        flat[k] = v['label']
                    else:
                        flat = flat | flatten_crystal_cdict(v)
                elif k == 'atoms' or k == 'simulation_cell' or k == 'job_details':
                    for i in v:
                        flat[i['label']] = i['value']
                elif k == 'software':
                    flat[k] = v[0]['label']
                else:
                    flat[k] = v
        return flat

def flatten_cdict(cdict):
        flat = {}
        for k, v in cdict.items():
            if k != '@context':
                if isinstance(v, dict):
                    if 'label' in v.keys():
                        flat[k] = v['label']
                    else:
                        flat = flat | flatten_cdict(v)
                elif isinstance(v, list):
                    for i in v: 
                        if isinstance(i, dict):
                            try:
                                flat[i['label']] = i['value']
                            except KeyError: # silently skips over terms that do not have label, value keys
                                pass
                        else:
                            flat[k] = v
                elif k == 'software':
                    flat[k] = v[0]['label']
                else:
                    flat[k] = v
        return flat

def GenericLammpsJobObject(o, space, project, collection, concept_dict, struct_concept_dict = None):
    # o = openbis_login(url, user)

    # if not space:
    #     space = user.upper()
    
    ob_coll = '/' + space + '/' + project + '/' + collection
    # cdict = flatten_lammps_cdict(concept_dict)
    cdict = flatten_cdict(concept_dict)

    objects = o.get_objects(
        space      = space,
        type       ='PYIRON_JOB.LAMMPS',
        start_with = 0,
    )
    exists = False
    for object_ in objects:
        if object_.p.get('$name') == cdict['job_name']:
            exists = True
            found_object = object_

    if exists:
        print("===================\n")
        print(f"Job already exists! Found job in: {found_object.identifier}\n")
        print("===================\n")
        print("Found job properties:\n")
        return found_object.p
    
    else:

        if cdict['job_status'] == 'finished':
            job_status = True
        else:
            job_status = False

        from datetime import datetime
        delta = datetime.strptime(cdict['job_stoptime'], "%Y-%m-%d %H:%M:%S") - datetime.strptime(cdict['job_starttime'], "%Y-%m-%d %H:%M:%S")

        atom_calc_type = None
        if 'molecular_statics' in concept_dict.keys() and 'minimization_algorithm' in cdict.keys():
            atom_calc_type = ('ATOM_CALC_STRUC_OPT').lower()
            description = 'Lammps simulation using pyiron for energy minimization/structural optimization.<p><span style="color:hsl(240,75%,60%);">' + \
                        '<strong>Scroll down below other properties to view conceptual dictionary with ontological ids of selected properties and values.</strong></span>' + \
                        '<br>The conceptual dictionary is in JSON-LD format. Learn more about it <a href="https://www.w3.org/ns/json-ld/">here</a></p>'
        
        atom_ionic_min_algo = None
        if cdict['minimization_algorithm'] == 'fire':
            atom_ionic_min_algo = ('MIN_ALGO_FIRE').lower()

        json_file = cdict['project_name'] + '/' + cdict['job_name'] + '_concept_dict.json'
        with open(json_file, 'r') as file:
            json_string = file.read()
        json_string = format_json_string(json_string)

        props_dict = {
            '$name': cdict['job_name'],
            'description': description,
            'workflow_manager': cdict['workflow_manager'],
            # 'bam_username': user,
            'bam_username': space.lower(),
            'sim_job_finished': job_status,
            'start_date': cdict['job_stoptime'],
            'sim_walltime_in_hours': delta.total_seconds()/3600,
            'sim_coretime_in_hours': cdict['sim_coretime_hours'],
            'ncores': cdict['number_cores'],
            'atomistic_calc_type': atom_calc_type,
            'periodic_boundary_x': cdict['periodicity_in_x'],
            'periodic_boundary_y': cdict['periodicity_in_y'],
            'periodic_boundary_z': cdict['periodicity_in_z'],
            'atom_cell_vol_relax': True if 'http://purls.helmholtz-metadaten.de/asmo/CellVolumeRelaxation' in cdict['dof'] else False,
            'atom_cell_shp_relax': True if 'CellShapeRelaxation' in cdict['dof'] else False,
            'atom_pos_relax': True if 'http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation' in cdict['dof'] else False,
            'atom_ionic_min_algo': atom_ionic_min_algo,
            'max_iters': cdict['maximum iterations'],
            'atom_e_tol_ion_in_ev': cdict['ionic energy tolerance'],
            'atom_f_tol_in_ev_a': cdict['force tolerance'],
            'atom_ionic_steps': cdict['NumberIonicSteps'],
            'atom_fin_tot_eng_in_ev': cdict['FinalTotalEnergy'],
            'atom_fin_vol_in_a3': cdict['FinalTotalVolume'],
            'atom_fin_pot_eng_in_ev': cdict['FinalPotentialEnergy'],
            'atom_force_max_in_ev_a': cdict['FinalMaximumForce'],
            'conceptual_dictionary': json_string,
        }
    
        object_ = o.new_object(
            type       = 'PYIRON_JOB.LAMMPS',
            space      = space,
            experiment = ob_coll,
            props      = props_dict
        )
        object_.save()

        #hdf_ver = job.to_dict()['HDF_VERSION']
        path_to_h5 = cdict['project_name'] + '/' + str(cdict['job_name']) + '.h5'
        path_to_json = cdict['project_name'] + '/' + str(cdict['job_name']) + '_concept_dict.json'
        path_to_yml = cdict['project_name'] + '/' + str(cdict['job_name']) + '_environment.yml'
        
        dataset_props_dict = {
            '$name': cdict['job_name'] + '.h5',
            'production_date': datetime.strptime(cdict['job_stoptime'], "%Y-%m-%d %H:%M:%S").date().strftime("%Y-%m-%d"),
            'file_format': 'HDF5',
            #'hdf5_version': hdf_ver,
            'reference': 'https://github.com/pyiron/pyiron_atomistics/blob/main/pyiron_atomistics/lammps/base.py',
        }
        
        ds_hdf = o.new_dataset(
            type       = 'PYIRON_JOB',
            collection = ob_coll,
            object     = object_,
            files      = [path_to_h5],
            kind = 'PHYSICAL',
            props      = dataset_props_dict
        )

        ds_hdf.save()

        ds_json = o.new_dataset(
            type       = 'ATTACHMENT',
            collection = ob_coll,
            object     = object_,
            files      = [path_to_json],
            kind = 'PHYSICAL',
            props      = {'$name':cdict['job_name']+'_concept_dict.json'}
        )

        ds_json.save()

        ds_yml = o.new_dataset(
            type       = 'COMP_ENV',
            collection = ob_coll,
            object     = object_,
            files      = [path_to_yml],
            kind = 'PHYSICAL',
            props      = {'$name':cdict['job_name']+'_environment.yml', 'env_tool': 'conda'}
        )

        ds_yml.save()

        if struct_concept_dict != None:
            # struct_object_name = 'input_structure_' + flatten_crystal_cdict(struct_concept_dict)['job_name']
            struct_object_name = 'input_structure_' + flatten_cdict(struct_concept_dict)['job_name']  ####################
            objects_1 = o.get_objects(
                space      = space,
                type       ='MAT_SIM_STRUCTURE.CRYSTAL',
                start_with = 0,
            )
            exists = False
            for object_1 in objects_1:
                if object_1.p.get('$name') == struct_object_name:
                    exists = True
                    found_object_1 = object_1

            if exists == True:
                object_.set_parents(found_object_1.identifier)
                object_.save()
            else:
                print("==============================\n")
                print("Create structure object first!\n")
                print("==============================\n")
        
        return object_.p

#TBD
def GenericCrystalObject(o, space, project, collection, concept_dict):
    # o = openbis_login(url, user)

    # if not space:
    #     space = user.upper()

    ob_coll = '/' + space + '/' + project + '/' + collection
    # cdict = flatten_crystal_cdict(concept_dict)
    cdict = flatten_cdict(concept_dict)

    ob_project_obj = o.get_project('/' + space + '/' + project)
    objects = ob_project_obj.get_objects(type = 'MAT_SIM_STRUCTURE.CRYSTAL')
    exists = False
    for object_ in objects:
        if object_.p.get('$name') == cdict['structure_name']:
            exists = True
            found_object = object_

    if exists == True:
        print("=========================\n")
        print(f"Structure already exists in project! Found structure in: {found_object.identifier}\n")
        print("=========================\n")
        print("Found structure properties:\n")
        return found_object.p

    else:

        description =   'Crystal structure generated using pyiron.<p><span style="color:hsl(240,75%,60%);">' + \
                        '<strong>Scroll down below other properties to view conceptual dictionary with ontological ids of selected properties and values.</strong></span>' + \
                        '<br>The conceptual dictionary is in JSON-LD format. Learn more about it <a href="https://www.w3.org/ns/json-ld/">here</a></p>'

        #description =   '<figure class="image image_resized image-style-align-left" style="width:9%;">' + \
        #                '<img src="/openbis/openbis/file-service/eln-lims/30/77/c8/3077c8f8-0e55-41e5-9021-7ed43863e5a4/cc44149e-b2d4-4850-9f88-7a4cd54c1a66.png">' + \
        #                '</figure><p><br>Lammps simulation using pyiron for energy minimization/structural optimization.<br>&nbsp;</p><p><span style="color:hsl(240,75%,60%);">' + \
        #                '<strong>Scroll down below other properties to view conceptual dictionary with ontological ids of selected properties and values.</strong></span>' + \
        #                '</p><p>The conceptual dictionary is in JSON-LD format. Learn more about it <a href="https://www.w3.org/ns/json-ld/">here</a><br>&nbsp;</p>'

        
        species = {}
        for i in concept_dict['atoms']:
            if i['label'] != 'total_number_atoms':
                species[i['label']] = i['value']

        json_file = cdict['path'] + cdict['structure_name'] + '_concept_dict.json'
        with open(json_file, 'r') as file:
            json_string = file.read()
        json_string = format_json_string(json_string)

        props_dict = {
            '$name': cdict['structure_name'],
            'description': description,
            'workflow_manager': cdict['workflow_manager'],
            'chem_species_by_n_atoms': str(species),
            'n_atoms_total': cdict['total_number_atoms'],
            'sim_cell_lengths_in_a': cdict['simulation_cell_lengths'],
            'sim_cell_vectors': cdict['simulation_cell_vectors'],
            'sim_cell_angles_in_deg': cdict['simulation_cell_angles'],
            'sim_cell_volume_in_a3': cdict['simulation_cell_volume'],
            'conceptual_dictionary': json_string,
        }

        if 'crystal_orientation' in cdict.keys(): props_dict['crystal_orientation'] = cdict['crystal_orientation']
        if 'lattice_parameter_a' in cdict.keys(): props_dict['lattice_param_a_in_a'] = cdict['lattice_parameter_a']
        if 'lattice_parameter_b' in cdict.keys(): props_dict['lattice_param_b_in_a'] = cdict['lattice_parameter_b']
        if 'lattice_parameter_c' in cdict.keys(): props_dict['lattice_param_c_in_a'] = cdict['lattice_parameter_c']
        if 'lattice_parameter_c_over_a' in cdict.keys(): props_dict['lattice_c_over_a'] = cdict['lattice_parameter_c_over_a']
        if 'lattice_angle_alpha' in cdict.keys(): props_dict['lattice_angalpha_in_deg'] = cdict['lattice_angle_alpha']
        if 'lattice_angle_beta' in cdict.keys(): props_dict['lattice_angbeta_in_deg'] = cdict['lattice_angle_beta']
        if 'lattice_angle_gamma' in cdict.keys(): props_dict['lattice_anggamma_in_deg'] = cdict['lattice_angle_gamma']
        if 'lattice_volume' in cdict.keys(): props_dict['lattice_volume_in_a3'] = cdict['lattice_volume']

        if 'space_group' in cdict.keys():
            spg_map = get_space_group_mapping(cdict['space_group'])
            props_dict['space_group'] = spg_map
        if 'bravais_lattice' in cdict.keys():
            bvl_map = get_bravais_lattice_mapping(cdict['bravais_lattice'])
            props_dict['bravais_lattice'] = bvl_map

        object_ = o.new_object(
            type       = 'MAT_SIM_STRUCTURE.CRYSTAL',
            space      = space,
            experiment = ob_coll,
            props      = props_dict
        )
        object_.save()

        struct_file_name = cdict['structure_name'] + '.h5'
        struct_file_path = cdict['path'] + cdict['structure_name'] + '.h5'
        upload_atomistic_structure_file(o, object_, struct_file_name, ob_coll, struct_file_path)
        cdict_file_name = cdict['structure_name'] + '_concept_dict.json'
        cdict_file_path = cdict['path'] + cdict['structure_name'] + '_concept_dict.json'
        upload_concept_dict(o, object_, cdict_file_name, ob_coll, cdict_file_path)
    
    return object_.p

def get_space_group_mapping(spg):
    if spg == 'Im-3m':
        return ('SPACE_GROUP.IM-3M').lower()
    elif spg == 'Fm-3m':
        return ('SPACE_GROUP.FM-3M').lower()
    elif spg == 'P6_3/mmc':
        return ('SPACE_GROUP.P63_MMC').lower()
    else:
        raise ValueError(f'Invalid Bravais lattice, maybe a formatting error?')

def get_bravais_lattice_mapping(bvl):
    if bvl == 'bcc':
        return ('BODY_CENTER_CUBIC').lower()
    elif bvl == 'fcc':
        return ('FACE_CENTER_CUBIC').lower()
    elif bvl == 'hcp':
        return ('HEX_CLOSE_PACK').lower()
    else:
        raise ValueError(f'Invalid Bravais lattice, maybe a formatting error?')

def upload_atomistic_structure_file(o, structure_object, structure_name, collection, file_path):
    dataset_props_dict = {
        '$name': structure_name,
        'multi_mat_scale': 'Electronic/Atomistic',
        'sw_compatibility': 'ASE',
        'file_format': 'HDF5',
    }
    
    ds_hdf = o.new_dataset(
        type       = 'MAT_SIM_STRUCTURE',
        collection = collection,
        object     = structure_object,
        files      = [file_path],
        kind       = 'PHYSICAL',
        props      = dataset_props_dict
    )
    ds_hdf.save()
    return

def upload_concept_dict(o, ob_object, cdict_name, collection, file_path):
    ds_json = o.new_dataset(
        type       = 'ATTACHMENT',
        collection = collection,
        object     = ob_object,
        files      = [file_path],
        kind       = 'PHYSICAL',
        props      = {'$name': cdict_name}
        )
    ds_json.save()
    return

