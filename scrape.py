import os
from os.path import isdir
from itertools import product
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-inpt', type=str) #input dir
parser.add_argument('-out', type=str) #output dir
parser.add_argument('-overwrite', default=True) #overwrite output files?
args = parser.parse_args()


dirs = [item for item in os.scandir(args.inpt) if isdir(item)]
keys = ["Active site"]

def scrape(fname,f):
    ''' scrape html file for subunits, active site residues
    '''
    # turn html into list of strings
    f = [str(line).split('>') for line in f]
    # ^^ gives [[$,$,$,...]], f[0] is [$,$,$,...]
#    print(f)
#3

    kln = []
    # attempt to add the enzyme name to the list containing the residue numbers of
    # the active site indices
    for i,line in enumerate(f):
        for ii,block in enumerate(line):
            if '<h1 property="name"' in block:
                #print('blockline',line[ii+1:ii+3])
                kln.append(line[ii+1].split('<')[0])
                break

    for i,line in enumerate(f):
        # each line is a list of blocks from a string split with '>''
        for ii,block in enumerate(line):
            # each block is a string of text seperated by
            for k in keys:
                #print(i,line)
                if k in block:
                    residue_number = line[ii+1].split('<')[0]


                    # residue numeber should be able to be an integer
                    try:
                        residue_number = int(residue_number)
                        kln.append(residue_number)
                    except:
                        pass
    return kln


def get_indices(info, mact):
    ''' info is a packet of the pdbid of the prot and
    the subunits it contains '''
    filename = info[0] + '.pdb' # ex: 1NEY.pdb
    units = info[1:]
    pdb = './{}/{}'.format(args.inpt,filename)

    def index(chain,oi,mact):
        act=[]

        with open (pdb,'r') as f:

            for ii,line in enumerate(f):
                line = line.split(' ')
                line = list(filter(lambda x: x!='',line))

                if (
                    line[0] == 'ATOM'
                    and line[4] == chain
                    and int(line[5]) in oi
                    ):
                        act.append(ii)

            mact.append(act)
        return mact


    def find_chain(substructure, subunit_number, oi, mact):
        # get the name of the chain from the pdb file

        ''''subunit' may be either term depending on how it
        appears in the html. pdb format capitalizes. This just
        capitalizes the str for searching pdb therefor.
        '''
        substructure = substructure.upper()

        combo = substructure + ' ' + subunit_number
        # ex: combo == 'CHAIN B'

        with open(pdb,'r') as f:

            # lines in pdb file containing word "COMPND"
            compound_lines = [l for l in f if "COMPND" in l]

            for i,line in enumerate(compound_lines):

                if combo in line:

                    if "MOLECULE" in line:
                        chain_line = compound_lines[i+1]
                    elif "SYNONYM" in line:
                        chain_line = compound_lines[i-1]

                    if not isinstance(chain_line, list):
                        chain_line = chain_line.split(' ')

                    for ii,c in enumerate(chain_line):
                        if c == 'CHAIN:':
                            chain = chain_line[ii+1]
                            # clean up punctuation
                            chain = chain.replace(';','')
                            chain = chain.replace(',','')
                            mact = index(chain,oi,mact)

    new = []

    # units are subunits
    # iterate over each subunit

    def search_through_multiple():

            if len(u) > 1:
                su = u[0].split(' ')
                original_indices = u[1:]

                for i,word in enumerate(su):

                    # subunit/chain number/letter follows these words in file
                    if (word == 'subunit' or word == 'chain'):
                        #print(word, u, info[0])
                        substructure = word

                        try:
                            subunit_number = su[i+1] # <-- see prev. comment
                            # get rid of empy spaces in str
                            subunit_number = subunit_number.replace(',','')

                            # grab the chain/subunit from the file
                            #try:
                            chain = find_chain(substructure,
                                            subunit_number,
                                            original_indices,
                                            mact
                                    )
                            #except FileNotFoundError as e:
                            #    print('FileNotFoundError', os.strerror(e.errno))
                            #    pass
                        except IndexError:
                            pass
    #print(pdb)


    #def search_through_one():

    act=[]
    oi = units[0][1:]
    #print(oi)

    with open (pdb,'r') as f:

        for i,line in enumerate(f):
            line = line.split(' ')
            line = list(filter(lambda x: x!='',line))

            try:

                if (
                    line[0] == 'ATOM'
                    and int(line[5]) in oi
                    ):
                        print(line[5], oi)
                        act.append(i)

            # it seems that one of the files has an atomic COORDINATE
            # at position 5 in the line
            # this the exception and passing
            except ValueError:
                pass
        mact.append(act)


    #search_through_one()

    return mact

#########################################################
####################  MAIN ##############################
#########################################################


def master():
    print('start')
    ''' This function exists solely to save space below
    as it would otherwise be written twice'''

    #list of lists of residue numbers, 1 sub-list per subunit
    protein_as_residues = [dir.name]
    # 'meta active site indices' -- 1 per protein
    mact = []

    for file in os.scandir(dir):
        # just gives name of PDB ID (omitting '.pdb' from file name)
        fname = file.name[:-4]

        try:

            with open(file, 'r') as f:
                # search pdb file for active site atom lines
                # subunit_as_residues are for each subunit
                subunit_as_residues = scrape(fname,f)

                # don't add empty lists for subunits w/o as indices
                if not subunit_as_residues == []:
                    protein_as_residues.append(subunit_as_residues)


        except IsADirectoryError:
            # sometimes dirs are found in input dir
            pass

    # if the list of AS residues contains more than just the header
    if len(protein_as_residues) > 1:

        mact = get_indices(protein_as_residues,mact)

        # don't bother writing to file if no indices found
        if not mact == []:

            """ print out list of line numbers of atoms
            belonging to residues in the active site"""
            for x in mact:
                cmd = 'echo {} >> {}'.format(x,ind_out)
                print(cmd)
                os.system(cmd)


for dir in dirs:

    ''' "try, except" is a means of setting the output dir name to input
    dir name by default unless args.out is given
    '''
    try:
        # if output dir doesn't exist
        if not os.path.isdir(args.out):
            # make it
            os.system('mkdir ./{}'.format(args.out))

        outpt_pdb_dir = './{}/{}/'.format(args.out,dir.name)

        if not os.path.isdir(outpt_pdb_dir):
            cmd = 'mkdir ' + outpt_pdb_dir
            #print(cmd)
            os.system(cmd)

        # init 1 dir and 1 file per protein to fill with as indices later
        ind_out = './{}/{}/{}.index'.format(args.out,dir.name,dir.name) # for outputting AS line index file later

    except:
        # will be caught if args.out doesn't exist, and default to args.inpt
        ind_out = './{}/{}/{}.index'.format(args.inpt,dir.name,dir.name)

    if args.overwrite is True:
        # old output files will be overwritten
        if os.path.isfile(ind_out):
            os.system('rm {}'.format(ind_out))
            os.system('touch {}'.format(ind_out))
            master()
        # if not file (is dir)
        else:
            master()
    elif args.overwrite is False:
        # old output files will not be overwritten
        if os.path.isfile(ind_out):
            pass
        else:
            master()
