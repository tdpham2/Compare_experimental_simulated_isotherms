import json
import sys
import numpy as np
import glob
import subprocess
import os
from ase.io import read as ase_read

def resample(pressures, uptakes):
    """ Resample list of pressures and uptakes. If less than 10 pressures per isotherm, then don't need to do anything
    Parameters:
        pressures: list of pressures, in Pa
        uptakes: uptakes
    Returns: 
        spressures: resampled pressures
        suptakes: resampled uptakes
    
    """
    def get_closest_value(v, list_input):
        """ Return the index of closest value to v from a list of input
        Parameters:
            v: float
            list_input: list of float
        """
        abs_diff = [abs(v-i) for i in list_input]
        abs_diff = np.array(abs_diff)

        ind = np.argmin(abs_diff)
        return ind

    Pmax = max(pressures)
    if len(pressures) != len(uptakes):
        print("Pressures and Uptakes have different length")
        return 'Error'
    if Pmax < 110000:
        sampled_pressures = [5000, 10000, 15000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000]
    if Pmax >= 110000:
        sample_dist = int((Pmax-10000)/12)
        sampled_pressures = list(range(10000, int(Pmax)+1, sample_dist))
    spressures = []
    suptakes = []

    for p in sampled_pressures:
        ind = get_closest_value(p, pressures)
        if pressures[ind] not in spressures:
            spressures.append(pressures[ind])
            suptakes.append(uptakes[ind])
    return spressures, suptakes

def convert_to_Pa(unit, pressures):
    """ Convert a list of pressures to Pa
    Parameters:
        unit: Unit of the input. Accept "bar", "kPa", "Pa", "atm"
        pressure: A list of float
    """
    if unit == 'bar':
        c_pressures = [round(i*100000) for i in pressures]
    if unit == 'kPa':
        c_pressures = [round(i*1000)   for i in pressures]
    if unit == 'Pa':
        c_pressures = [round(i) for i in pressures]
    if unit == 'atm':
        c_pressures = [round(i*101325) for i in pressures]

    return c_pressures

def get_isotherm_json(fname):
    with open(fname, 'r') as f:
        data = json.load(f)

    adsorbate_name = data['adsorbates'][0]['name']
    adsorbent_name = data['adsorbent']['name']
    units = data['adsorptionUnits']
    isotherm_data = data['isotherm_data']
    pressure_unit = data['pressureUnits']
    temperature = data['temperature']
    pressures = [item['pressure'] for item in isotherm_data]
    uptakes = [item['species_data'][0]['adsorption'] for item in isotherm_data]
    adsorption_units = data['adsorptionUnits']

    return adsorbate_name, temperature, pressure_unit, pressures, uptakes, adsorption_units, adsorbent_name

def create_raspa_input(cif, output_dir, pressure, adsorbate, temperature=298, fname='simulation.input', charges=True, cycles=10000, cutoff=12.8, raspa_input_path=None):
    """ Create RASPA2 simulation.input file
    Parameters: 
        cif: str, path to cif file
        output_dir: str, directory to write simulation.input file
        pressure: list, list of pressure in Pa
        adsorbate: str, adsorbate name, best to use the same name as the .def file
        temperature: float, temperature in K. Default is 298 K
        fname: str, input file name. Default is simulation.input
        charges: bool, specify whether to use charges. Default: True
        cycles: int, number of initialization and equilibration cycles.
        cutoff: float, VDW and Coulomb cutoff.
    """
    def calculate_cell_size(cif, cutoff):
        cif = open(cif,"r")
        mof_cif= cif.read()
        cif.close()

        for line in mof_cif.split("\n"):
            if "_cell_length_a" in line:
                length_a = line.split()[1]  #unit cell vector
                length_a =float(length_a)  #string to float
            if "_cell_length_b" in line:
                length_b = line.split()[1]
                length_b = float(length_b)
            if "_cell_length_c" in line:
                length_c= line.split()[1]
                length_c= float(length_c)
            if "_cell_angle_alpha" in line:
                alpha = line.split()[1]
                alpha = float(alpha)
            if "_cell_angle_beta" in line:
                beta= line.split()[1]
                beta= float(beta)
            if "_cell_angle_gamma" in line:
                gamma = line.split()[1]
                gamma = float(gamma)
        ax = length_a
        ay = 0.0
        az = 0.0
        bx = length_b * np.cos(gamma * np.pi / 180.0)
        by = length_b * np.sin(gamma * np.pi / 180.0)
        bz = 0.0
        cx = length_c * np.cos(beta * np.pi / 180.0)
        cy = (length_c * length_b * np.cos(alpha * np.pi /180.0) - bx * cx) / by
        cz = (length_c ** 2 - cx ** 2 - cy ** 2) ** 0.5
        unit_cell =  np.asarray([[ax, ay, az],[bx, by, bz], [cx, cy, cz]])
        #Unit cell vectors
        A = unit_cell[0]
        B = unit_cell[1]
        C = unit_cell[2]
        #minimum distances between unit cell faces
        Wa = np.divide(np.linalg.norm(np.dot(np.cross(B,C),A)), np.linalg.norm(np.cross(B,C)))
        Wb = np.divide(np.linalg.norm(np.dot(np.cross(C,A),B)), np.linalg.norm(np.cross(C,A)))
        Wc = np.divide(np.linalg.norm(np.dot(np.cross(A,B),C)), np.linalg.norm(np.cross(A,B)))

        uc_x = int(np.ceil(cutoff/(0.5*Wa)))
        uc_y = int(np.ceil(cutoff/(0.5*Wb)))
        uc_z = int(np.ceil(cutoff/(0.5*Wc)))

        return uc_x, uc_y, uc_z

    uc_x, uc_y, uc_z = calculate_cell_size(cif, cutoff)

    adsorbates = {"hydrogen":"H2", "carbon dioxide":"CO2", "nitrogen":"N2", "water": "H2O", "argon":"Ar", "xenon": "Xe", "krypton": "Kr", "methane":"CH4"}

    if '/' in cif:
        cifname = cif.split('/')[-1][:-4]
    else:
        cifname = cif[:-4]
    if len(pressure) == 1:

        with open(output_dir + '/' + fname, 'w') as f:
            f.write("SimulationType                MonteCarlo\n")
            f.write("NumberOfCycles                {}\n".format(cycles))
            f.write("NumberOfInitializationCycles  {}\n".format(cycles))
            f.write("PrintEvery                    1000\n")
            f.write("PrintForcefieldToOutput       no\n")
            f.write("Movies                        Yes\n")
            f.write("WriteMoviesEvery              {}\n".format(int(cycles/20)))
            f.write("Cutoff                        {}\n".format(cutoff))
            f.write("Forcefield                    Local\n")
            f.write("ChargeMethod                  Ewald\n")
            f.write("EwaldPrecision                1e-6\n")
            f.write("CutoffCoulomb                 {}\n".format(cutoff))

            if charges == True:
                f.write("UseChargesFromCIFFile         yes\n")
            else:
                f.write("UseChargesFromCIFFile         no\n")

            f.write("Framework                     0\n")
            f.write("FrameworkName                 {}\n".format(cifname))
            f.write("UnitCells                     {} {} {}\n".format(uc_x, uc_y, uc_z))
            f.write("ExternalTemperature           {}\n".format(temperature))
            f.write("ExternalPressure              {}\n".format(pressure[0]))

            f.write("Component 0  MoleculeName             {}\n".format(adsorbate))
            f.write("             MoleculeDefinition       Local\n")
            f.write("             TranslationProbability   0.5\n")
            f.write("             RotationProbability      0.5\n")
            f.write("             ReinsertionProbability   0.5\n")
            f.write("             SwapProbability          1.0\n")
            f.write("             CreateNumberOfMolecules  0\n")

    else:
        for p in pressure:
            subprocess.run("mkdir -p {}/{}".format(output_dir, str(p)), shell=True)
            subprocess.run("cp {} {}/{}".format(cif, output_dir, str(p)), shell=True)
            with open(output_dir + '/' + str(p) + '/' + fname, 'w') as f:
                f.write("SimulationType                MonteCarlo\n")
                f.write("NumberOfCycles                {}\n".format(cycles))
                f.write("NumberOfInitializationCycles  {}\n".format(cycles))
                f.write("PrintEvery                    1000\n")
                f.write("PrintForcefieldToOutput       no\n")
                f.write("Movies                        Yes\n")
                f.write("WriteMoviesEvery              {}\n".format(int(cycles/20)))
                f.write("Cutoff                        {}\n".format(cutoff))
                f.write("Forcefield                    Local\n")
                f.write("ChargeMethod                  Ewald\n")
                f.write("EwaldPrecision                1e-6\n")
                f.write("CutoffCoulomb                 {}\n".format(cutoff))

                if charges == True:
                    f.write("UseChargesFromCIFFile         yes\n")
                else:
                    f.write("UseChargesFromCIFFile         no\n")

                f.write("Framework                     0\n")
                f.write("FrameworkName                 {}\n".format(cifname))
                f.write("UnitCells                     {} {} {}\n".format(uc_x, uc_y, uc_z))
                f.write("ExternalTemperature           {}\n".format(temperature))
                f.write("ExternalPressure              {}\n".format(p))

                f.write("Component 0  MoleculeName             {}\n".format(adsorbate))
                f.write("             MoleculeDefinition       Local\n")
                f.write("             TranslationProbability   0.5\n")
                f.write("             RotationProbability      0.5\n")
                f.write("             ReinsertionProbability   0.5\n")
                f.write("             SwapProbability          1.0\n")
                f.write("             CreateNumberOfMolecules  0\n")

            if raspa_input_path != None:
                if adsorbate == 'H2':
                    if temperature < 200:
                        adsorbate_path = raspa_input_path + '/' + 'H2_FH'
                    else:
                        adsorbate_path = raspa_input_path + '/' + 'H2'
                else:
                    adsorbate_path = raspa_input_path + '/' + adsorbate

                if os.path.isdir(adsorbate_path):
                    subprocess.run("cp {}/* {}/{}".format(adsorbate_path, output_dir, str(p)), shell=True)

dirs = sorted(glob.glob('/home/tdpham/work/research/NIST_ISODB/csd_calculations/subset1/H2/*/'))
raspa_input = '/home/tdpham/software/raspa_inputs'
adsorbates = {"Hydrogen":"H2", "Carbon Dioxide":"CO2", "Nitrogen":"N2", "Methane": "CH4"}

#isotherms = ['10.1039C1cc12858b.isotherm5.json', '10.1039C1cc12858b.isotherm6.json']
isotherms = ['10.1039C2cc34482c.Isotherm1.json']
path_to_isotherms = '/home/tdpham/software/Compare_experimental_simulated_isotherms/csd/matched_CIFs_and_DOI'
path_to_simulations = '/home/tdpham/work/research/Experimental_isotherm_project/simulated_isotherm'
raspa_input_path = 'RASPA_input'

cif_assignment = []
iso_assignment = []

with open('cif_assignment.txt', 'r') as f:
    for l in f:
        line = l.strip().split(',')
        iso_assignment.append(line[0])
        cif_assignment.append(line[1])

for iso in iso_assignment:
    doi = iso.lower().split('isotherm')[0][:-1].lower()
    path_to_iso = os.path.join(path_to_isotherms, doi, iso)
    if not os.path.isfile(path_to_iso):
        print("Error setting up {}".format(iso))
    else: 
        adsorbate_name, temperature, pressure_unit, pressures, uptakes, adsorption_units, adsorbent_name = get_isotherm_json(path_to_iso)
        adsorbate_label = adsorbates[adsorbate_name]
        output_dir = os.path.join(path_to_simulations, doi, iso[:-5])
        if os.path.isdir(output_dir):
            continue
        else:
            print(doi)
            cif_index = iso_assignment.index(iso)
            cif = os.path.join('CIFs', cif_assignment[cif_index])
            if os.path.isfile(cif):
                converted_pressures = convert_to_Pa(pressure_unit, pressures)
                spressures, suptake = resample(converted_pressures, uptakes)
                if adsorbate_label == 'CH4':
                    create_raspa_input(cif, output_dir, spressures, adsorbates[adsorbate_name], temperature=temperature, fname='simulation.input', charges=False, cycles=10000, cutoff=12.8, raspa_input_path=raspa_input_path)

                else:
                    create_raspa_input(cif, output_dir, spressures, adsorbates[adsorbate_name], temperature=temperature, fname='simulation.input', charges=True, cycles=10000, cutoff=12.8, raspa_input_path=raspa_input_path)
            else:
                print("Error setting up {}. No CIF found!!!".format(iso))
"""
for d in dirs:
    print(d)
    jsons = glob.glob(d+'/'+'*.json')
    cifs = glob.glob(d+'/'+'*.cif')
    js_name = [i.split('/')[-1] for i in jsons]
    js_dir = [i[:-5] for i in jsons]
    if len(cifs) == 1:
        cif = cifs[0]
    # TODO: Create cif assignment for directories with multiple CIFs
    else:
        js_assign = []
        cif_assign = []

        if os.path.isfile(d + 'cif_assignment.txt'):
            with open(d + 'cif_assignment.txt', 'r') as f:
                for line in f:
                    data = line.strip().split(',')
                    js_assign.append(data[0])
                    cif_assign.append(data[1])
        else:
            #print("Missing cif_assigment.txt at {}".format(d))
            continue
    for index, i in enumerate(js_name):
        if len(cifs) > 1:
            try:
                ind_assign = js_assign.index(i)
            except ValueError:
                continue
            cif_path = glob.glob(d + '/' + cif_assign[ind_assign] + '*.cif')
            if len(cif_path) > 1:
                print('Error with {}'.format(d))
                break
            else:
                cif = cif_path[0]
        #subprocess.run("mkdir -p {}".format(i), shell=True)
        #subprocess.run("cp {} {}".format(cif, i), shell=True)

        adsorbate_name, temperature, pressure_unit, pressures, uptakes, adsorption_units, adsorbent_name = get_isotherm_json(jsons[index])
        adsorbate_name = adsorbate_name.lower()
        converted_pressures = convert_to_Pa(pressure_unit, pressures)
        spressures, suptake = resample(converted_pressures, uptakes)
        create_raspa_input(cif, js_dir[index], spressures, adsorbates[adsorbate_name], ndir=len(converted_pressures), temperature=temperature, fname='simulation.input', charges=False, cycles=10000, cutoff=14.0, raspa_input_path=raspa_input)

"""
