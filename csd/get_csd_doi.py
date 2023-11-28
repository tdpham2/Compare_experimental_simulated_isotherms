from ccdc import io

csd_ref = []

with open('csd_mof.csv','r') as f:
    for line in f:
        csd_ref.append(line.strip())

for index, ref in enumerate(csd_ref):
    csd_reader = io.EntryReader('CSD')
    try:
        crys = csd_reader.crystal(ref)
        entry = csd_reader.entry(ref)
        doi = entry.publication.doi.lower()
    except RuntimeError:
        with open('CSD_missing_REF.txt', 'a') as f:
            f.write('{}\n'.format(ref))
            continue
    except AttributeError:
        with open('CSD_missing_DOI.txt', 'a') as f:
            f.write('{}\n'.format(ref))
            doi = None
    with open('csd_with_doi.csv','a') as f:
        f.write("{},{}\n".format(ref, doi))
