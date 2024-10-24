"""
LAMMPS specific functions for parsing

Use this a reference for specific implementations
"""
import os
import numpy as np
import ast
import json


def process_lammps_job(job):
    method_dict = {}
    add_contexts(method_dict)
    #get_structures(job, method_dict)
    identify_method(job, method_dict)
    extract_calculated_quantities(job, method_dict)
    add_software(job, method_dict)
    #add_pyiron_details(job, method_dict)
    #add_sim_details(job, method_dict)
    get_simulation_folder(job, method_dict)

    file_name = job.path + '_concept_dict.json'
    with open(file_name, 'w') as f:
        json.dump(method_dict, f, indent=2)
    return method_dict

def process_structure_crystal(job, initial_structure = None, rotation_indices = None):
    sample_dict = {}
    add_contexts_sample(sample_dict)
    get_chemical_species(job, sample_dict)
    get_simulation_cell(job, sample_dict)
    add_software(job, sample_dict, struct=True)
    get_simulation_folder(job, sample_dict)

    struct_file_name = job.path + '_input_structure_concept_dict.json'
    with open(struct_file_name, 'w') as f:
        json.dump(sample_dict, f, indent=2)
    return sample_dict

def add_contexts(method_dict):
    method_dict['@context'] = {}
    method_dict['@context']['sample'] = 'http://purls.helmholtz-metadaten.de/cmso/AtomicScaleSample'
    method_dict['@context']['path'] = 'http://purls.helmholtz-metadaten.de/cmso/hasPath'
    method_dict['@context']['dof'] = 'http://purls.helmholtz-metadaten.de/asmo/hasRelaxationDOF'
    method_dict['@context']['inputs'] = 'http://purls.helmholtz-metadaten.de/asmo/hasInputParameter'
    method_dict['@context']['label'] = 'http://www.w3.org/2000/01/rdf-schema#label'
    method_dict['@context']['unit'] = 'http://purls.helmholtz-metadaten.de/asmo/hasUnit'
    method_dict['@context']['value'] = 'http://purls.helmholtz-metadaten.de/asmo/hasValue'
    method_dict['@context']['outputs'] = 'http://purls.helmholtz-metadaten.de/cmso/hasCalculatedProperty'
    #method_dict['@context']['workflow_manager'] = ''
    #method_dict['@context']['software'] = ''
    method_dict['@context']['molecular_dynamics'] = 'http://purls.helmholtz-metadaten.de/asmo/MolecularDynamics'
    method_dict['@context']['molecular_statics'] = 'http://purls.helmholtz-metadaten.de/asmo/MolecularStatics'
    method_dict['@context']['ensemble'] = 'http://purls.helmholtz-metadaten.de/asmo/hasStatisticalEnsemble'

def add_contexts_sample(sample_dict):
    sample_dict['@context'] = {}
    sample_dict['@context']['path'] = 'http://purls.helmholtz-metadaten.de/cmso/hasPath'
    sample_dict['@context']['unit_cell'] = 'http://purls.helmholtz-metadaten.de/cmso/UnitCell'
    sample_dict['@context']['atoms'] = 'http://purls.helmholtz-metadaten.de/cmso/Atom'
    sample_dict['@context']['molecules'] = 'http://purls.helmholtz-metadaten.de/cmso/Molecule'
    sample_dict['@context']['bravais_lattice'] = 'http://purls.helmholtz-metadaten.de/cmso/hasBravaisLattice'
    sample_dict['@context']['chemical_species'] = 'http://purls.helmholtz-metadaten.de/cmso/ChemicalSpecies'
    sample_dict['@context']['simulation_cell'] = 'http://www.w3.org/2000/01/rdf-schema#label'
    sample_dict['@context']['label'] = 'http://www.w3.org/2000/01/rdf-schema#label'
    sample_dict['@context']['unit'] = 'http://purls.helmholtz-metadaten.de/cmso/hasUnit'
    sample_dict['@context']['value'] = 'http://purls.helmholtz-metadaten.de/asmo/hasValue'
    sample_dict['@context']['vector'] = 'http://purls.helmholtz-metadaten.de/cmso/Vector'
    #method_dict['@context']['workflow_manager'] = ''
    #method_dict['@context']['software'] = ''

def get_simulation_folder(job, method_dict):
    method_dict['path'] = os.path.join(job.project.path, f'{job.name}_hdf5')

def get_structures(job, method_dict):
    initial_pyiron_structure = job.structure
    final_pyiron_structure = job.get_structure(frame=-1)
    
    method_dict['sample'] =  {'initial':initial_pyiron_structure, 
                            'final': final_pyiron_structure}

def identify_method(job, method_dict):
    job_dict = job.input.to_dict()
    input_dict = {
        job_dict["control_inp/data_dict"]["Parameter"][x]: job_dict[
            "control_inp/data_dict"
        ]["Value"][x]
        for x in range(len(job_dict["control_inp/data_dict"]["Parameter"]))
    }
    dof = []
    temp = None
    press = None
    md_method = None
    ensemble = None

    if "min_style" in input_dict.keys():
        dof.append("http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation")
        if job.input.control['fix___ensemble'] != 'all nve': dof.append("http://purls.helmholtz-metadaten.de/asmo/CellVolumeRelaxation")
        md_method = "molecular_statics"
        e_tol = float(input_dict['minimize'].split()[0])
        f_tol = float(input_dict['minimize'].split()[1])
        maxiter = int(input_dict['minimize'].split()[2])

    elif "nve" in input_dict["fix___ensemble"]:
        if int(input_dict["run"]) == 0:
            md_method = "molecular_statics"
            ensemble = "http://purls.helmholtz-metadaten.de/asmo/MicrocanonicalEnsemble"

        elif int(input_dict["run"]) > 0:
            dof.append("http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation")
            md_method = "molecular_dynamics"
            ensemble = "http://purls.helmholtz-metadaten.de/asmo/MicrocanonicalEnsemble"

    elif "nvt" in input_dict["fix___ensemble"]:
        raw = input_dict["fix___ensemble"].split()
        temp = float(raw[3])
        dof.append("http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation")
        md_method = "molecular_dynamics"
        ensemble = "http://purls.helmholtz-metadaten.de/asmo/CanonicalEnsemble"

    elif "npt" in input_dict["fix___ensemble"]:
        dof.append("http://purls.helmholtz-metadaten.de/asmo/AtomicPositionRelaxation")
        dof.append("http://purls.helmholtz-metadaten.de/asmo/CellVolumeRelaxation")
        if "aniso" in input_dict["fix___ensemble"]:
            dof.append("http://purls.helmholtz-metadaten.de/asmo/CellShapeRelaxation")
        md_method = "molecular_dynamics"
        raw = input_dict["fix___ensemble"].split()
        temp = float(raw[3])
        press = float(raw[7])
        ensemble = "http://purls.helmholtz-metadaten.de/asmo/IsothermalIsobaricEnsemble"

    method_dict[md_method] = {}

    if md_method == "molecular_statics":
        method_dict[md_method]['minimization_algorithm'] = input_dict['min_style']

    input_dict = {
        job_dict["control_inp/data_dict"]["Parameter"][x]: job_dict[
            "control_inp/data_dict"
        ]["Value"][x]
        for x in range(len(job_dict["control_inp/data_dict"]["Parameter"]))
                }
    pb = []
    pb.append(input_dict['boundary'])
    if (pb[0][0] == "p"):
        method_dict[md_method]['periodicity_in_x'] = True
    else:
        method_dict[md_method]['periodicity_in_x'] = False
    if (pb[0][2] == "p"):
        method_dict[md_method]['periodicity_in_y'] = True
    else:
        method_dict[md_method]['periodicity_in_y'] = False
    if (pb[0][4] == "p"):
        method_dict[md_method]['periodicity_in_z'] = True
    else:
        method_dict[md_method]['periodicity_in_z'] = False

    method_dict[md_method]['inputs'] = []

    temperature = {}
    temperature["value"] = temp
    temperature["unit"] = "K"
    temperature["label"] = "temperature"

    method_dict[md_method]['inputs'].append(temperature)

    pressure = {}
    pressure["value"] = press
    pressure["unit"] = "GigaPA"
    pressure["label"] = "pressure"

    method_dict[md_method]['inputs'].append(pressure)

    energy_tol = {}
    energy_tol["value"] = e_tol
    energy_tol["unit"] = "EV"
    energy_tol["label"] = "ionic energy tolerance"

    method_dict[md_method]['inputs'].append(energy_tol)

    force_tol = {}
    force_tol["value"] = f_tol
    force_tol["unit"] = "EV-PER-ANGSTROM"
    force_tol["label"] = "force tolerance"

    method_dict[md_method]['inputs'].append(force_tol)

    maximum_iterations = {}
    maximum_iterations["value"] = maxiter
    maximum_iterations["label"] = "maximum iterations"

    method_dict[md_method]['inputs'].append(maximum_iterations)

    method_dict[md_method]["ensemble"] = ensemble
    
    method_dict["dof"] = dof


    # now process potential
    inpdict = job.input.to_dict()
    ps = inpdict["potential_inp/data_dict"]["Value"][0]
    name = inpdict["potential_inp/potential/Name"]
    potstr = job.input.to_dict()["potential_inp/potential/Citations"]
    potdict = ast.literal_eval(potstr[1:-1])
    url = None
    if "url" in potdict[list(potdict.keys())[0]].keys():
        url = potdict[list(potdict.keys())[0]]["url"]

    if 'meam' in ps:
        method_dict['@context']['potential'] = "http://purls.helmholtz-metadaten.de/asmo/ModifiedEmbeddedAtomModel"
    elif 'eam' in ps:
        method_dict['@context']['potential'] = "http://purls.helmholtz-metadaten.de/asmo/EmbeddedAtomModel"
    elif 'lj' in ps:
        method_dict['@context']['potential'] = "http://purls.helmholtz-metadaten.de/asmo/LennardJonesPotential"
    elif 'ace' in ps:
        method_dict['@context']['potential'] = "http://purls.helmholtz-metadaten.de/asmo/MachineLearningPotential"
    else:
        method_dict['@context']['potential'] = "http://purls.helmholtz-metadaten.de/asmo/InteratomicPotential"


    method_dict[md_method]["potential"] = {}
    method_dict[md_method]["potential"]["label"] = name
    if url is not None:
        method_dict[md_method]["potential"]["@id"] = url

def add_software(job, method_dict, struct=False):
    method_dict["workflow_manager"] = {}
    method_dict["workflow_manager"]["@id"] = "http://demo.fiz-karlsruhe.de/matwerk/E457491"
    import subprocess
    import platform
    try:
        if "Windows" in platform.system():
            output1 = subprocess.check_output(['findstr', 'pyiron_atomistics', job.project.name + '\\' + job.name + '_environment.yml'])
        else:
            output1 = subprocess.check_output(['grep', 'pyiron_atomistics', job.project.name + '/' + job.name + '_environment.yml'])
        s1 = str((output1.decode('utf-8')))
    except:
        s1 = ''
    try:
        if "Windows" in platform.system():
            output2 = subprocess.check_output(['findstr', 'pyiron_workflow', job.project.name + '\\' + job.name + '_environment.yml'])
        else:
            output2 = subprocess.check_output(['grep', 'pyiron_workflow', job.project.name + '/' + job.name + '_environment.yml'])
        s2 = str((output2.decode('utf-8')))
    except:
        s2 = ''
    # hdf_ver = job.to_dict()['HDF_VERSION']
    st = 'p' + s1.split('=')[0].split('p')[1] + "=" + s1.split('=')[1] + ', ' + 'p' + s2.split('=')[0].split('p')[1] + "=" + s2.split('=')[1] #+ ', pyiron_HDF_version=' + hdf_ver
    method_dict["workflow_manager"]["label"] = st

    pyiron_job_details = []
    pyiron_job_details.append(
        {
            "label": "job_name",
            "value": job.name,
        }
    )
    pyiron_job_details.append(
        {
            "label": "project_name",
            "value": job.project.name,
        }
    )
    if not struct:
        pyiron_job_details.append(
            {
                "label": "job_type",
                "value": job.database_entry.hamilton,
            }
        )
        pyiron_job_details.append(
            {
                "label": "job_status",
                "value": str(job.status),
            }
        )
        pyiron_job_details.append(
            {
                "label": "job_starttime",
                "value": str(job.database_entry.timestart.strftime("%Y-%m-%d %H:%M:%S")),
            }
        )
        pyiron_job_details.append(
            {
                "label": "job_stoptime",
                "value": str(job.database_entry.timestop.strftime("%Y-%m-%d %H:%M:%S")),
            }
        )
        pyiron_job_details.append(
            {
                "label": "sim_coretime_hours",
                "value": job.database_entry.totalcputime,
            }
        )
        pyiron_job_details.append(
            {
                "label": "number_cores",
                "value": job.to_dict()['server']['cores'],
            }
        )
 
        software = {
            "@id": "http://demo.fiz-karlsruhe.de/matwerk/E447986",
            "label": "LAMMPS " + job.to_dict()['executable']['version'],
        }
        method_dict["software"] = [software]

    method_dict["job_details"] = pyiron_job_details

def extract_calculated_quantities(job, method_dict):
    """
    Extracts calculated quantities from a job.

    Parameters
    ----------
    job : pyiron.Job
        The job object containing the calculated quantities.

    Returns
    -------
    list
        A list of dictionaries, each containing the label, value, unit, and associate_to_sample of a calculated quantity.

    """
    aen = np.mean(job.output.energy_tot)
    fen = job.output.energy_tot[-1]
    fpe = job.output.energy_tot[-1]
    avol = np.mean(job.output.volume)
    fvol = job.output.volume[-1]
    fmax = job.output.force_max[-1]
    nionic = len(job.output.steps)-1
    outputs = []
    outputs.append(
        {
            "label": "AverageTotalEnergy",
            "value": np.round(aen, decimals=4),
            "unit": "EV",
        }
    )
    outputs.append(
        {
            "label": "FinalTotalEnergy",
            "value": np.round(fen, decimals=4),
            "unit": "EV",
        }
    )
    outputs.append(
        {
            "label": "FinalPotentialEnergy",
            "value": np.round(fpe, decimals=4),
            "unit": "EV",
        }
    )
    outputs.append(
        {
            "label": "AverageTotalVolume",
            "value": np.round(avol, decimals=4),
            "unit": "ANGSTROM3",
        }
    )
    outputs.append(
        {
            "label": "FinalTotalVolume",
            "value": np.round(fvol, decimals=4),
            "unit": "ANGSTROM3",
        }
    )
    outputs.append(
        {
            "label": "FinalMaximumForce",
            "value": np.round(fmax, decimals=16),
            "unit": "EV-PER-ANGSTROM",
        }
    )
    outputs.append(
        {
            "label": "NumberIonicSteps",
            "value": nionic,
        }
    )
    method_dict['outputs'] =  outputs

def identify_unit_cell(job, sample_dict):
    #Stuff lattice quantities and orientation
    pass

def identify_crystal_structure(job, sample_dict):
    #Stuff here for space group and bravais lattice
    pass

def get_chemical_species(job, sample_dict):
    structure = job.structure
    natoms = structure.get_number_of_atoms()
    species_dict = dict(structure.get_number_species_atoms())
    atoms_list = []
    for k in species_dict.keys():
        element = {}
        element["value"] = species_dict[k]
        element["label"] = k
        atoms_list.append(element)

    atoms_list.append({'value': natoms, 'label': 'total_number_atoms'})
        
    sample_dict['atoms'] = atoms_list

def get_simulation_cell(job, sample_dict):
    structure = job.structure
    cell_lengths = str([structure.cell.cellpar()[0],structure.cell.cellpar()[1],structure.cell.cellpar()[2]])
    cell_vectors = str([structure.cell[0],structure.cell[1],structure.cell[2]])
    cell_angles = str([structure.cell.cellpar()[3],structure.cell.cellpar()[4],structure.cell.cellpar()[5]])
    cell_volume = structure.get_volume()
    
    simulation_cell_details = []

    cell_lengths_dict = {}
    cell_lengths_dict["value"] = cell_lengths
    cell_lengths_dict["unit"] = "ANGSTROM"
    cell_lengths_dict["label"] = 'simulation_cell_lengths'
    cell_lengths_dict["@id"] = "http://purls.helmholtz-metadaten.de/cmso/hasLength"
    simulation_cell_details.append(cell_lengths_dict)
    
    cell_vector_dict = {}
    cell_vector_dict["value"] = cell_vectors
    cell_vector_dict["unit"] = "ANGSTROM"
    cell_vector_dict["label"] = 'simulation_cell_vectors'
    cell_vector_dict["@id"] = "http://purls.helmholtz-metadaten.de/cmso/hasVector"
    simulation_cell_details.append(cell_vector_dict)

    cell_angles_dict = {}
    cell_angles_dict["value"] = cell_angles
    cell_angles_dict["unit"] = "DEGREES"
    cell_angles_dict["label"] = 'simulation_cell_angles'
    cell_angles_dict["@id"] = "http://purls.helmholtz-metadaten.de/cmso/hasAngle"
    simulation_cell_details.append(cell_angles_dict)

    cell_volume_dict = {}
    cell_volume_dict["value"] = np.round(cell_volume, decimals=4)
    cell_volume_dict["unit"] = "ANGSTROM3"
    cell_volume_dict["label"] = 'simulation_cell_volume'
    cell_volume_dict["@id"] = "http://purls.helmholtz-metadaten.de/cmso/hasVolume"
    simulation_cell_details.append(cell_volume_dict)

    sample_dict['simulation_cell'] = simulation_cell_details

