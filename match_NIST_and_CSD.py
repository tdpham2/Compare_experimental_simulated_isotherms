from ccdc import io
import sys
import pandas as pd
import subprocess

nist_doi = pd.read_csv("isodb-library/DOI_mapping.csv", sep=', ', skiprows=0, names=['doi_normal','dir'])
csd_doi = pd.read_csv("csd/csd_with_doi.csv", sep=',', names=['ref', 'doi'])
csd_doi.dropna(inplace=True)

csd_outpath = 'csd/matched_CIFs_and_DOI' 
nist_iso = 'isodb-library/Library/'
nist_doi['doi'] = nist_doi['doi_normal'].str.lower()

nist_dir_list = nist_doi['dir'].tolist()
nist_doi_list = nist_doi['doi'].tolist()
csd_doi_list = csd_doi['doi'].tolist()
csd_ref_list = csd_doi['ref'].tolist()

csd_reader = io.EntryReader('CSD')

for count, doi in enumerate(nist_doi_list):
    print(doi)
    match = [index for index, i in enumerate(csd_doi_list) if i == doi]
    refs = [csd_ref_list[i] for i in match]
    if len(refs) == 0:
        print(doi)
        continue
    else:
        subprocess.run('cp -r {}/{} {}'.format(nist_iso, nist_dir_list[count], csd_outpath), shell=True)
        for ref in refs:
            crys = csd_reader.crystal(ref)
            with io.CrystalWriter("{}/{}/{}.cif".format(csd_outpath, nist_dir_list[count], ref)) as mol_writer:
                mol_writer.write(crys)
