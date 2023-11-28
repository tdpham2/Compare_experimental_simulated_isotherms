import json
import sys
import glob
import subprocess
import os

all_adsorbates = {'Methane': 'CH4', 'Hydrogen': 'H2', 'Nitrogen': 'N2', 'Oxygen': 'O2', 'Carbon Dioxide': 'CO2', 'Benzene': 'Benzene', 'Carbon monoxide': 'CO', 'Water': 'H2O', 'Ethane': 'C2H6', 'Ethene': 'C2H4', 'Methanol': 'CH3OH', 'Argon': 'Ar', '2-Propanol': '2-Propanol', '1-Butanol': '1-Butanol', 'Ethanol': 'C2H5OH', 'Acetone': 'Acetone', '1-Propanol': '1-Propanol', 'Propene': 'Propene', 'N-propane': 'N-propane', 'Acetylene': 'Acetylene', 'Deuterium': 'D2', 'Xenon': 'Xe', 'Krypton': 'Kr', 'Acetonitrile': 'Acetonitrile', 'Nitric oxide': 'NO', '1-Pentanol': '1-Pentanol', 'N-Hexane': 'N-Hexane', 'Cyclohexane': 'Cyclohexane', 'P-Xylene': 'P-Xylene', 'N-Butane': 'N-Butane', '2-Methylbutane': '2-Methylbutane', 'N-Pentane': 'N-Pentane', 'Isobutane': 'Isobutane', 'Dimethyl ether': 'Dimethyl-ether', 'Toluene': 'Toluene', 'Isobutanol': 'Isobutanol', 'Pyridine': 'Pyridine', 'Tert-Butanol': 'Tert-Butanol', 'Helium': 'Helium', 'Mercury Dichloride': 'Mercury-Dichloride', 'Hydrogen sulfide': 'H2S', 'O-Xylene': 'O-Xylene', 'M-Xylene': 'M-Xylene', 'Ammonia': 'NH3', 'Ethylbenzene': 'Ethylbenzene', 'Trichloromethane': 'Trichloromethane', 'Diethyl ether': 'Diethyl-ether', 'Neon': 'Ne', 'Nitrobenzene': 'Nitrobenzene', 'Sulfur dioxide': 'SO2', 'Tetrahydrofuran': 'Tetrahydrofuran', 'Acetaldehyde': 'Acetaldehyde'}

nist_path = 'csd/matched_CIFs_and_DOI'
output = 'simulated_isotherms/'

nist_iso = glob.glob(nist_path + '/*/')

for item in nist_iso:
    cifs = glob.glob(item + '/*.cif')
    cifnames = [i.split('/')[-1].split('.')[0] for i in cifs]

    isotherms = glob.glob(item + '/' + '*.json')
    cifs = glob.glob(item + '/' + '*.cif')
    doi = item.split('/')[-2]
    dict_adsorbents = {}
    dict_adsorbates = {}
    for iso in isotherms:
        f = open(iso)
        data = json.load(f)
        f.close()
        adsorbates = data['adsorbates']
        if len(adsorbates) == 1:
            adsorbate_name = data['adsorbates'][0]['name']
        else:
            print(iso)
            continue
        subprocess.run('mkdir -p {}/{}'.format(output, all_adsorbates[adsorbate_name]), shell=True)
        subprocess.run('mkdir -p {}/{}/{}'.format(output, all_adsorbates[adsorbate_name], doi), shell=True)
        subprocess.run('cp {} {}/{}/{}'.format(iso, output, all_adsorbates[adsorbate_name], doi), shell=True)
        for c in cifs:
            subprocess.run('cp {} {}/{}/{}'.format(c, output, all_adsorbates[adsorbate_name], doi), shell=True)
    
