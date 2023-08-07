from Bio.PDB import *
import pymol2
import os

def remove_duplicates(pdb, pdb_code, parser):
    ''' Remove disordered atoms -- whose positions are alternatives to
    primary positions in pdb file. Save new pdb with no duplicate atoms
    '''

    class NotDisordered(Select):
        def accept_atom(self, atom):
            return atom.get_altloc() in (' ','A')

    # get structure
    struct = parser.get_structure("my_pdb", pdb)

    # avoid occupancy warning from BioPython by setting occ
    for atom in struct.get_atoms():
        atom.set_occupancy(1)

    # save tailored structure
    io = PDBIO()
    io.set_structure(struct)
    file_wo_duplicates = pdb_code + '_no_dup.pdb'
    io.save(file_wo_duplicates, select=NotDisordered())

    # return file name
    return file_wo_duplicates


def add_hydrogen(pdb, pdb_code):
    ''' Add hydrogen to pdb file.
    '''
    pdb_with_hydrogen = pdb_code + '_wH.pdb'
    with pymol2.PyMOL() as pymol:
        pymol.cmd.load(pdb, pdb_code)
        pymol.cmd.h_add()
        pymol.cmd.save(pdb_with_hydrogen)
    return pdb_with_hydrogen


def preprocess_pdb(pdb, pdb_code, parser):
    ''' Prepare PDB for analysis.
    '''
    cleaned_pdb = remove_duplicates(pdb, pdb_code, parser)
    processed_pdb = add_hydrogen(cleaned_pdb, pdb_code)

    # remove temp file
    os.system('rm {}'.format(cleaned_pdb))

    # return new file from which model will be read
    return processed_pdb


def get_solvated_pdb(pdb):
    ''' Get path of solvated pdb file if called for by
    args.solvated in main py file.
    '''
    split_path = os.path.split(pdb)
    pdb_code, ext = split_path[-1].split('.')
    solvated_filename = pdb_code + '_solvated.' + ext

    # ideally, pdb is located in input_pdbs/
    if len(split_path) > 1:
        path_to = '/'.join(split_path[:-1]) + '/'
        return path_to + solvated_filename
    # ...but may be located in main dir
    else:
        return solvated_filename
