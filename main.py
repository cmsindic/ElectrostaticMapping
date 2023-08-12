import os
from argparse import ArgumentParser
from Residues import *
import ElectricField as ef
from fileProcessing.PDB import get_solvated_pdb
from spatial.cartesian import magnitude
from time import time

parser = ArgumentParser(
    prog='Electrostatic Field Collector',
    description="""Map electrostatic fields in a pdb file
    from amber99sb forcefield without needing gromacs or other
    such software to ascertain charges from forcefield files""")

parser.add_argument('--filename', required=True)
parser.add_argument('--cutoff', type=float, required=False)
parser.add_argument('--save', action='store_true', required=False)
parser.add_argument('--indexdir', type=str,
                    required=False, default='index_files')
parser.add_argument('--activesite', action='store_true', required=False)
parser.add_argument('--grid', action='store_true',required=False)
parser.add_argument('--solvated', action='store_true',required=False)
args = parser.parse_args()

original_pdb = args.filename
pdb_code = os.path.split(original_pdb)[-1][:-4]
index_file_dir = args.indexdir
index_file = os.path.join(index_file_dir, pdb_code + '.index')
rtp_file = 'amber99sb-ildn.ff/aminoacids.rtp'

''' Get residues from solvated pdb if pdb is solvated using PACKMOL
or some other method, else get residues from original pdb.'''
pdb = original_pdb if not args.solvated else get_solvated_pdb(original_pdb)
model = ResiduesInFile(pdb, pdb_code, rtp_file)
model.get_residues()

# (charge, array([x, y, z])) for all atoms in pdb file
cac = model.get_charges_and_coords(model.residues)
charges, coords = tuple(zip(*cac))

''' Need to collect active site information from original file
given that the line numbers of the active atoms are taken from
the file before solvation.'''
active_site = ActiveSite(index_file, original_pdb)
active_site.get_residues(residues=model.residues)
active_site.get_coords()

t1 = time()

model.convert_solvent_to_water_objects()

in_enzyme = set(
    w for w in model.solvent
        if target_near_coords(
            target=w.oxygen.coord,
            coords=model.solute_coords,
            cutoff=4)
)

near_as = set(
    w for w in model.solvent
        if target_near_coords(
            target=w.oxygen.coord,
            coords=active_site.coords,
            cutoff=4)
)

bulk_solvent = set(model.solvent) - in_enzyme
in_enzyme -= near_as

def get_field_across_sele(sele, charges_coords):
    sele_h = [h for water in sele for h in water.hydrogen]
    return ef.get_avg_e_field(targets=sele_h,
                              charges_coords=charges_coords,
                              exclude=[],
                              cutoff=10)

partitions = {'Bulk Solvent    ': bulk_solvent,
              'Near Enzyme     ': in_enzyme,
              'Near Active Site': near_as}

print("Average electric field strength on water protons:")
print("Sample Partition  |  Field Strength (a.u.)  |  N Hydrogen")
print("---------------------------------------------------------")
for k, v in partitions.items():
    sele = list(v)
    avg_magnitude = magnitude(
        get_field_across_sele(sele, cac)
        )
    avg_magnitude = round(avg_magnitude, 7)
    print(k, " | ",
          avg_magnitude,
          "             | ",
          len(sele)
        )


print("Run time:", time()-t1)




'''
for water in model.solvent:
    for atom in water.hydrogen:

        field = ef.electric_field(
            target=atom.coord,
            charges_coords=local_coords,
            exclude=water.coords)
'''

'''
RES = 0.5
if args.grid is True:
    from spatial.cartesian import *
    from time import time

    t1 = time()

    xyz = []
    for cube in generate_cubes(coords, res=RES, padding=1):
        cube.check_if_contains_point(coords)
        if cube.contains_point:
            xyz.append(cube.get_central_point())


    from matplotlib import pyplot as plt
    x, y, z = list(zip(*xyz))
    plt.scatter(x, y)
    plt.show()

    print(time()-t1)
'''
