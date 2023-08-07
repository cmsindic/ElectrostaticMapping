import os
import scandir
from sys import argv

base = 'https://www.uniprot.org/uniprot/'


def fill(fpath, fname):

    print(fpath, fname)
    prot = fname[:-4]
    url = base + prot
    cmd = ('wget --output-document={} {}').format(fpath, url)

    print("Running // " + cmd)
    os.system(cmd)


for item in scandir.scandir(argv[1]):
    if os.path.isdir(item.path):
        for file in scandir.scandir(item.path):
            os.system('rm {} && touch {}'.format(file.path,file.path))
            fill(file.path, file.name)
