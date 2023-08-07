from spatial.cartesian import *

class Residue():
    ''' Extension of BIO PDB residue class created to rename hydrogen in PDB
    to match those is RTP and ultimately correlate atom charges with thier
    coordinates in the PDB file.
    '''
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
        ''' Get atomic partial charge for a given atom.
        '''
        for atom in self.atoms:
            for atom_name_in_rtp, __, charge, __ in self.atom_lines:
                if atom.name == atom_name_in_rtp:
                    atom.charge = charge
                    break
            else:
                atom.charge = None

    def get_bond_partners(self):
        ''' Get bond partners for a given atom.
        '''

        def get_partner(atom, bond_pair):
            ''' Get the other atom in the rtp bond pair line.
            '''
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

        def closest_atom_to(h):
            ''' Get the closest atom to a given hydrogen.
            '''

            def has_room(a):
                ''' Does the atom have room for more bonds to H?
                '''
                return a.n_children - len(a.bonded_h) > 0

            # available atoms
            av_atoms = [a for a in self.res_non_h if has_room(a)]
            # sort by distance to h
            av_atoms.sort(key=lambda x: dist(x.coord, h.coord))
            return av_atoms[0]

        # for each H in PDB res, find closest atom, assign bond
        # do until number of H assigned == number of H in RTP
        res_h = self.res_h
        n_h_to_be_bonded = len(self.rtp_hydrogen)
        for _ in range(n_h_to_be_bonded):
            for h in res_h:
                closest_atom = closest_atom_to(h)
                closest_atom.bonded_h.append(h)
                res_h.remove(h)
                break

        # replace name of H in file with name of H in bond
        for atom in self.res_non_h:
            for i,h in enumerate(atom.bonded_h):
                h.name = atom.rtp_hydrogen[i]
