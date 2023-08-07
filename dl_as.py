import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
import os
from time import sleep
from sys import argv
import scandir

# run like "python dl_as.py isomerases"
doi = './{}'.format(argv[1]) # dir of interest
generic_rcsb_url = "https://www.rcsb.org/structure/"
generic_uniprot_url = "https://www.uniprot.org/uniprot/"

print('Starting...')

def rm_ext(f):
    return (f[0][:-4], f[1])

files = [(f.name, f.path) for f in scandir.scandir(doi)]
files = list(filter(lambda x: '.pdb' in x[0], files))
files = list(map(lambda x: rm_ext(x), files))

def get_rcsb_html(f):
    fname, fpath = f
    rcsb_url = generic_rcsb_url + fname
    return rcsb_url

def get_uniprot_urls(url):
    blocks=[]
    r = requests.get(url)
    r = r.text
    html = list(r.split('<'))

    for block in html:
        if '.org/uniprot/' in block:
            blocks.append(block[-6:])

    for i,b in enumerate(blocks):
        blocks[i] = generic_uniprot_url + b

    return blocks

def get_indices(urls, fname):
    for url in urls:
        if not os.path.isdir(argv[1]):
            os.system("mkdir {}".format(argv[1]))
        if not os.path.isdir("{}/{}".format(argv[1],fname)):
            os.system("mkdir {}/{}".format(argv[1],fname))
        os.system("touch {}/{}/{}.txt".format(argv[1],fname, url[-6:]))

for file in files:
    print(file)

    to_write = '{}/{}'.format(argv[1],file[0])

    if not os.path.isdir(to_write):
        os.system('mkdir {}'.format(to_write))
    else:
        pass

    rcsb_url = get_rcsb_html(file)
    uni_urls = get_uniprot_urls(rcsb_url)

    #sleep(1)
    as_indices = get_indices(uni_urls, file[0])
    #sleep(1)
'''
base = 'https://www.uniprot.org/uniprot/'

def fill(fpath, fname):

    print(fpath, fname)
    prot = fname[:-4]
    url = base + prot
    curl = "$(curl --ciphers DEFAULT@SECLEVEL=1 {})"
    bcmd = "echo " + curl + " >> {}"
    cmd = bcmd.format(url, fpath)

    print("Running // " + cmd)
    os.system(cmd)


for item in scandir.scandir(argv[1]):
    if os.path.isdir(item):
        for file in scandir.scandir(item):
            fill(file.path, file.name)
            '''
