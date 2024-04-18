import os
import pymol2


# Collect PDB files in a list
model_dirs = [d for d in os.scandir() if os.path.isdir(d.name)]
input_files = [f.path for d in model_dirs for f in os.scandir(d)]


def add_hydrogen(pdb_path):
    ''' Add hydrogen to pdb file.
    Overwrite file with result.
    '''
    
    # Skip files already proccessed
    nametag = '_wH.pdb'
    if nametag in pdb_path:
        return
        
    print(pdb_path)

    # path/to/1234.pdb --> path/to/1234
    pdb_no_ext = os.path.splitext(pdb_path)[0]
    
    # path/to/1234 --> path/to/1234_wH.pdb
    pdb_with_hydrogen = pdb_no_ext + nametag
        
    pymol.cmd.load(pdb_path)
    pymol.cmd.h_add()
    pymol.cmd.save(pdb_with_hydrogen)
    
    return pdb_with_hydrogen


# Add Hydrogen to each PDB and save
for f in input_files:
    with pymol2.PyMOL() as pymol:
        try:
            add_hydrogen(f)
        except pymol.CmdException:
            pass