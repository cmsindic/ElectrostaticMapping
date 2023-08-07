import os
import Bio
from Bio.PDB import *
import pymol2
import numpy as np
from itertools import combinations

pdb = "1NEY.pdb"
pdbcode = pdb[:-4]


###################### GET INFORMATION FROM RTP FILE ######################


def fetch_rtp(rtp_file):
    '''Get lines with charges, atom names, etc.,
    from rtp_atoms.rtp. Return dicts with resnames
    as keys
    '''

    class RTPline():
        ''' To contain information about lines in RTP file.
        '''
        def __init__(self,line):
            self.is_resline = line[0] == '['
            self.is_bond_declaration = 'bond' in line
            self.is_dihedral = 'dihedrals' in line
            self.is_improper = 'impropers' in line
            self.is_start_of_atoms_section = 'atoms' in line
            self.is_end_of_bond_section = self.is_dihedral or self.is_improper
            line = line.replace('\t',' ')
            line = line.replace('\n',' ')
            line = line.split(' ')
            self.content = [x for x in line if x not in ('','[',']')]
            self.is_empty = len(self.content) == 0

    # Collect atom information
    atom_lines, bond_lines = {}, {}

    for i,line in enumerate(rtp_file):
        # Line with resname has no whitespace as first char
        # Starts with [
        if line[0] != '[':
            continue

        # Create key for residue
        resname = line.split(' ')[1]
        atom_lines[resname] = []
        bond_lines[resname] = []

        # Start saving lines to list to be dict value
        j = 0
        line_type = None # None until line section type is declared
        while i + j + 1 < len(rtp_file):
            j += 1

            # line to class
            line = RTPline(rtp_file[i + j])

            if line.is_empty:
                continue

            # break if hit resline or start of bonds section
            if line.is_end_of_bond_section or line.is_resline:
                break
            # Don't add line declaring atoms, go to next line
            elif line.is_start_of_atoms_section:
                line_type = 'ATOM'
                continue
            # Don't add line declaring bonds, go to next line
            elif line.is_bond_declaration:
                line_type = 'BOND'
                continue

            if line_type == 'ATOM':
                id, atom, charge, indx = line.content
                formatted_line = (id, atom, float(charge), int(indx))
                atom_lines[resname].append(formatted_line)
            elif line_type == 'BOND':
                bond_lines[resname].append(line.content)

            def add_formatted_line_to_dict(line,d,line_type):
                ...

        # Skip all lines added to dicts to look for new res start line
        i += j

    return atom_lines, bond_lines


# Generate Dicts of atoms and bonds from rtp
with open('amber99sb-ildn.ff/aminoacids.rtp','r') as f:
    rtp_atoms, rtp_bonds = fetch_rtp([line for line in f])


############################################################################

###################### PREPROCESS PDB FILE #################################

''' Remove disordered atoms -- whose positions are alternatives to
primary positions in pdb file. Save new pdb with no duplicate atoms
'''
parser = PDBParser(QUIET=True)
s = parser.get_structure("my_pdb", pdb)
io = PDBIO()
io.set_structure(s)

class NotDisordered(Select):
    def accept_atom(self, atom):
        return not atom.is_disordered() or atom.get_altloc() == "A"

file_wo_duplicates = pdbcode + '_no_dup.pdb'
io.save(file_wo_duplicates, select=NotDisordered())

''' Add hydrogen to pdb file. '''
processed_pdb_file = pdbcode + '_wH.pdb'
with pymol2.PyMOL() as pymol:
    pymol.cmd.load(file_wo_duplicates, 'myprotein')
    pymol.cmd.h_add()
    pymol.cmd.save(processed_pdb_file)


############################################################################

############ MATCH ATOMS IN RTP WITH PDB, ASSIGN CHARGES ###################


def dist(x, y):
    ''' Euclidian distance between cartesian coordinates. '''
    d = 0
    for i in range(3):
        d += (x[i] - y[i])**2
    return d**0.5

class Residue():

    def __init__(self, r):
        self.bio_res_class = r
        self.name = r.get_resname()
        self.id = r.id[1]
        self.atoms = list(r.get_atoms())
        self.res_h = [a for a in self.atoms if a.element == 'H']
        self.res_non_h = [a for a in self.atoms if a.element != 'H']
        self.is_not_water = self.name not in ('HOH', 'SOL')

    def get_dict_values(self, rtp_atoms_dict, rtp_bonds_dict):
        ''' Try to find the residue in the dicts of atoms and bonds
        created from lines the the rtp file. Identify if PDB residue
        is in both dictionaries.
        '''
        try:
            self.atom_lines = rtp_atoms_dict[self.name]
            self.bond_lines = rtp_bonds_dict[self.name]
            self.in_dict = True
        except:
            self.in_dict = False

    def is_good_residue(self):
        ''' Return whether residue is non-water residue that can be
        found in the rtp file (atom/bond dict).
        '''
        return self.is_not_water and self.in_dict

    def get_atom_charges(self):
        ''' Get atomic partial charge for a given atom. '''
        for atom in self.atoms:
            for atom_name_in_rtp,__,charge,__ in self.atom_lines:
                if atom.name == atom_name_in_rtp:
                    atom.charge = charge
                    break
            else:
                atom.charge = None

    def get_bond_partners(self):
        ''' Get bond partners for a given atom. '''

        def get_partner(atom, bond_pair):
            ''' Get the other atom in the rtp bond pair line. '''
            bond_pair = set(bond_pair)
            bond_pair.discard(atom.name)
            bonded_atom = list(bond_pair)[0]
            return bonded_atom

        # obtain names of atoms bonded to each atom, specifically hydrogen
        for atom in self.atoms:
            pairs = [line for line in self.bond_lines if atom.name in line]
            atom.bond_partners = [get_partner(atom, p) for p in pairs]
            atom.rtp_hydrogen = [a for a in atom.bond_partners if a[0]=='H']

    def sort_hydrogen(self):
        ''' Translate H names in pdb with H charges in
        bonds list from RTP to assign charges to pdb H.
        '''
        # get list of all H in residue, per RTP
        self.rtp_hydrogen = []
        for atom in self.res_non_h:
            self.rtp_hydrogen += atom.rtp_hydrogen
            atom.n_children = len(atom.rtp_hydrogen) # n H bonded to atom
            atom.bonded_h = [] # init empty array of H assigned to atom

        def has_room(a):
            ''' Does the atom have room for more bonds to H? '''
            return a.n_children - len(a.bonded_h) > 0

        # for each H in PDB res, find closest atom, assign bond
        # do until number of H assigned == number of H in RTP
        res_h = self.res_h
        n_h_to_be_bonded = len(self.rtp_hydrogen)
        for _ in range(n_h_to_be_bonded):
            for h in res_h:
                # available atoms
                av_atoms = [a for a in self.res_non_h if has_room(a)]
                av_atoms.sort(key=lambda x: dist(x.coord, h.coord))
                closest_atom = av_atoms[0]
                closest_atom.bonded_h.append(h)
                res_h.remove(h)
                break

        # replace name of H in file with name of H in bond
        for atom in self.res_non_h:
            for i,h in enumerate(atom.bonded_h):
                h.name = atom.rtp_hydrogen[i]

# reminder:
# processed_pdb_file is PDB with H added and alt positions removed
struct = parser.get_structure(pdbcode, processed_pdb_file)

struct_res = struct.get_residues()
all_residues = [Residue(r) for r in struct_res]
non_water_residues = [r for r in all_residues if r.is_not_water]

for residue in non_water_residues:
    if residue.name == 'HIS': residue.name = 'HIE'
    residue.get_dict_values(rtp_atoms, rtp_bonds)

residues = [r for r in non_water_residues if r.is_good_residue()]
bad_residues = [r.name for r in non_water_residues if not r in residues]

charge_and_coords = []
for residue in residues:
    # get H bonded to non-H
    residue.get_bond_partners()
    # rename H in pdb to H in RTP
    residue.sort_hydrogen()
    # get bond partners, now with new H names
    residue.get_bond_partners()
    # get charges from RTP, now that H have proper names
    residue.get_atom_charges()

    # generate list of coords and charges
    for atom in residue.atoms:
        if atom.is_disordered():
            ''' If charge is none and atom as bond partners,
            something is wrong '''
            if atom.charge is None:
                print(residue.name)
                assert len(atom.bond_partners) == 0
        else:
            charge_and_coords.append((atom.charge, atom.coord))

charge_and_coords = [c for c in charge_and_coords if c[0] != None]
#charge_and_coords.sort(key=lambda x: x[0], reverse=True)


############################################################################

############## GET ACTIVE SITE LOCATION FROM ORIGINAL PDB ##################

from parse import *

# get AS line numbers from index file
index_file = pdbcode + '.index'
with open(index_file,'r') as f:
    # there is only 1 line in index files
    file = [line for line in f]
    assert len(file) == 1
    active_site_line_numbers = str2list(file[0])

with open(pdb,'r') as f:
    pdb_lines = [line for line in f]

raw_as_lines = [pdb_lines[i] for i in active_site_line_numbers]
active_site_lines = [pdbline(x) for x in raw_as_lines]
as_residues_from_pdb = {(a[3], int(a[5])) for a in active_site_lines}
as_residues = [r for r in residues if (r.name, r.id) in as_residues_from_pdb]

as_coords = []
for r in as_residues:
    for a in r.atoms:
        as_coords.append(a.coord)

central_as_point = np.mean(as_coords, axis=0)

############################################################################

############## GENERATE SPHERE OF POINTS AROUND ACTIVE SITE ################


def create_sphere(coord, r, resolution=7):
    '''
    create sphere with center (cx, cy, cz) and radius r
    '''
    cx, cy, cz = coord
    phi = np.linspace(0, 2*np.pi, 2*resolution)
    theta = np.linspace(0, np.pi, resolution)
    theta, phi = np.meshgrid(theta, phi)
    r_xy = r*np.sin(theta)
    x = cx + np.cos(phi) * r_xy
    y = cy + np.sin(phi) * r_xy
    z = cz + r * np.cos(theta)
    return x, y, z

radii = [x/2 for x in range(2,22,1)]
spheres = {}
for rad in radii:
    x,y,z = create_sphere(central_as_point, rad)
    sphere_coords = []
    for i in range(len(x)):
        for j in range(len(x[i])):
            sphere_coords.append((x[i][j], y[i][j], z[i][j]))
    spheres[rad] = sphere_coords


############################################################################

############## GET E FIELD AT POINTS AROUND ACTIVE SITE ####################

# coulomb constant
C = 9 * 10**9
# coulombs per elementary charge
e = 1.602 * 10**-19
angstromsToMeters = 10**-10
vmToAu = 5.14 * 10**11

def magnitude(x):
    # magnitude in 3 dimensions of a len=3 array of distances
    return sum(x**2)**0.5

def electricField(target, charge_and_coords):
    ''' Get E field for each frame of simulation
    '''
    field = np.zeros((1,3))
    for i,c in enumerate(charge_and_coords):
        charge, coord = charge_and_coords[i]
        xyz = angstromsToMeters * (np.array(target) - coord)
        d = magnitude(xyz)
        E = C * e * (charge / d**2)
        EComponents = (E / d) * xyz / vmToAu
        field += EComponents
    return field

for r, sphere_coords in spheres.items():
    for s in sphere_coords:
        for _, coords in charge_and_coords:
            if dist(s, coords) < 2:
                sphere_coords.remove(s)
                break
        else:
            continue
    spheres[r] = sphere_coords

avg_mags = []
for r, sphere_coords in spheres.items():
    sum_sphere_field_mag = 0
    for s in sphere_coords:
        field = electricField(s, charge_and_coords)
        sum_sphere_field_mag += magnitude(field[0])
    avg_sphere_field_mag = sum_sphere_field_mag / len(sphere_coords)
    avg_mags.append(avg_sphere_field_mag)
    print(r, avg_sphere_field_mag)

with open(pdbcode + '.csv', 'w') as f:
    f.write('Plot of E field mag vs distance from AS \n')
    f.write('r(A), E Field Mag (a.u.)\n')
    for i,r in enumerate(spheres.keys()):
        f.write(str(r) + ',' + str(avg_mags[i]) + '\n')

'''
from matplotlib import pyplot as plt
plot = plt.figure()
plt.scatter(spheres.keys(),avg_mags)
plt.show()
'''
