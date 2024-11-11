def classic_structure(pr, structure, structure_name_prefix):
    structure_path = pr.name + '/'
    structure_name = structure_name_prefix + '_input_structure'
    
    from pyiron_base.storage.hdfio import FileHDFio
    hdf = FileHDFio(structure_path + structure_name + '.h5')
    structure.to_hdf(hdf)

    # Cannot guarantee this will be before manipulations hence skipping get_unit_cell_parameters
    from ob.concept_dict import process_structure_crystal
    struct_cdict = process_structure_crystal(pr, structure, structure_name, structure_path)
    
    return struct_cdict

def classic_lammps(lammps_job):
    from ob.concept_dict import process_lammps_job, export_env

    export_env(lammps_job.path)

    lammps_cdict = process_lammps_job(lammps_job)
    return lammps_cdict

def classic_murn(murn_job):
    from ob.concept_dict import export_env
    export_env(murn_job.path)

    from ob.concept_dict import process_murnaghan_job, process_lammps_job
    child_jobs_cdict = []
    for jobs in murn_job.iter_jobs():
        # from pyiron_base.storage.hdfio import FileHDFio
        import platform
        if "Windows" in platform.system():
            import shutil
            shutil.copy(murn_job.path + '_environment.yml', jobs.path + '_environment.yml')
        else:
            import os
            os.system('cp ' + murn_job.path + '_environment.yml ' + jobs.path + '_environment.yml')
        child_cdict = process_lammps_job(jobs)
        child_jobs_cdict.append(child_cdict)

    job_cdict = process_murnaghan_job(murn_job)

    return job_cdict, child_jobs_cdict

def classic_murn_equil_structure(murn_job):

    from pyiron_base.storage.hdfio import FileHDFio
    from ob.concept_dict import process_structure_crystal, get_unit_cell_parameters

    equil_structure = murn_job.get_structure()
    structure_name = murn_job.name + '_equilibrium_structure'
    structure_path = murn_job.project.name + '/'
    hdf_equil = FileHDFio(structure_path + structure_name + '.h5')
    equil_structure.to_hdf(hdf_equil)

    # Must be conventional unit cell or primitive cell
    # Provide the structure before repetition or any other manipulations
    structure_parameters = get_unit_cell_parameters(equil_structure)
    
    struct_cdict = process_structure_crystal(murn_job.project, equil_structure, structure_name, structure_path, structure_parameters)
        
    return struct_cdict

def upload_classic_pyiron(pr, job, o, space, project, collection):

    structure = job.structure
    if structure:
        
        from ob.openbis import GenericLammpsJobObject, GenericCrystalObject, MurnaghanJobObject

        # Project env file
        from ob.concept_dict import export_env
        export_env(pr.path + pr.name)

        struct_dict = classic_structure(pr, structure, structure_name_prefix=job.name)
        ob_structure = GenericCrystalObject(o, space, project, collection, struct_dict, show_object=False)

        if 'lammps' in job.to_dict()['TYPE']:
            job_cdict = classic_lammps(job)
            ob_job = GenericLammpsJobObject(o, space, project, collection, job_cdict, show_object=False)
        elif 'murn' in job.to_dict()['TYPE']:
            job_cdict, child_jobs_cdict = classic_murn(job)
            ob_job = MurnaghanJobObject(o, space, project, collection, job_cdict, child_jobs_cdict, show_object=False)

            equil_struct_dict = classic_murn_equil_structure(job)
            ob_equil_structure = GenericCrystalObject(o, space, project, collection, equil_struct_dict, show_object=False)

            ob_job.add_children([ob_equil_structure])
        else:
            job_type = job.to_dict()['TYPE']
            print(f'The {job_type} job type is not implemented for OpenBIS upload yet.')
            return
        
        ob_job.add_parents([ob_structure])
        ob_job.save()
        
    else:
        print('This job does not contain a structure and will not be uploaded. \
                Please add structure before trying to upload.')