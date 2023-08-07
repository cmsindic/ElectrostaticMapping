from Bio.PDB import *
import fileProcessing.RTP
import fileProcessing.parse as parse
from fileProcessing.PDB import preprocess_pdb
from spatial.cartesian import *
import warnings


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

    def get_dict_values(self, rtp_atoms_dict, rtp_bonds_dict, name, forgiving=False):
        ''' Try to find the residue in the dicts of atoms and bonds
        created from lines the the rtp file. Identify if PDB residue
        is in both dictionaries.
        '''
        assert rtp_atoms_dict.keys() == rtp_bonds_dict.keys()
        
        self.has_name_in_dict = False
        self.atom_nums_match_dict = False
        
        if name in rtp_atoms_dict.keys():
            '''If there aren't an equal number of atoms in
            the rtp file as in the model from the pdb,then
            the atom in the rtp is misnamed and is not
            technically in the dict. If 'forgiving' argument
            is given as True, ignore this requirement.
            '''
            rtp_entry = rtp_atoms_dict[self.name]
            same_atom_count = len(self.atoms) == len(rtp_entry)
            if any((same_atom_count, forgiving==True)):
                self.atom_lines = rtp_atoms_dict[self.name]
                self.bond_lines = rtp_bonds_dict[self.name]
                self.in_dict = True
                return             
                
        self.in_dict = False
        return

    def is_good_residue(self):
        ''' Return whether residue is non-water residue that can be
        found in the rtp file (atom/bond dict).
        '''
        return self.is_not_water and self.in_dict

    def get_bond_partners(self):
        ''' Get bond partners for a given atom.
        '''

        def rename_water_atoms(atom):
            ''' Water hydrogen need to be renamed, and can be done
            so predictably.
            '''
            if atom.element == 'O':
                atom.name = 'OW'
            elif atom.element == 'H':
                if '1' in atom.name:
                    atom.name = 'HW1'
                elif '2' in atom.name:
                    atom.name = 'HW2'
            return atom

        def get_partner(atom, bond_pair):
            ''' Get the other atom in the rtp bond pair line.
            '''
            bond_pair = set(bond_pair)
            bond_pair.discard(atom.name)
            bonded_atom = list(bond_pair)[0]
            return bonded_atom

        if not self.is_not_water:
            self.atoms = [rename_water_atoms(a) for a in self.atoms]

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

    def get_atom_charges(self):
        ''' Get atomic partial charge for a given atom.
        '''
        
        
        def assign_charge(atom, op, rtp_atoms):
            ''' Assign atom charge, rename atom, and
            and update rtp_atoms by remove the atom upon
            matching charges/names.
            '''
            for line in rtp_atoms:
                atom_name_in_rtp, _, charge, _ = line
                # Try seeing if model atom has the
                # same name as the rtp atom
                if op == 'EQ':
                    cond = atom_name_in_rtp == atom.name
                # For atoms that haven't exact names in rtp,
                # check if part of the atom name is in any of 
                # the rtp atom names or vice versa
                elif op == "CONTAINS":
                    c1 = atom_name_in_rtp in atom.name
                    c2 = atom.name in atom_name_in_rtp
                    cond = c1 or c2
                if cond:
                    atom.name = atom_name_in_rtp
                    atom.charge = charge
                    rtp_atoms.remove(line)
                    return atom, rtp_atoms
            return atom, rtp_atoms

        rtp_atoms = self.atom_lines.copy()
        for atom in self.atoms:
            atom.charge = None
            atom, rtp_atoms = assign_charge(atom, 'EQ', rtp_atoms)
        
        for atom in self.atoms:
            if atom.charge is None:
                atom, rtp_atoms = assign_charge(atom, 'CONTAINS', rtp_atoms)
        
        self.get_bond_partners()
        self.sort_hydrogen()
        
        for atom in self.atoms:
            if atom.charge is None:
                atom, rtp_atoms = assign_charge(atom, 'CONTAINS', rtp_atoms)
        
        no_charges = [a.name for a in self.atoms if a.charge == None]
        assigned_atoms = [a.name for a in self.atoms if a.charge != None]
        '''
        if len(no_charges) > 0:
            print(self.name,len(self.atom_lines),len(self.atoms))
            print(no_charges)
            print(rtp_atoms)'''

        '''
        if len(no_charges) != 0:
            print(self.name, self.id)
            print(no_charges)
            print([a[0] for a in self.atom_lines if not a[0] in assigned_atoms])
            print()
        '''




class Water(Residue):
    ''' Special methods for handling water if electrostatic fields
    are to be measured as they are felt by water's components.
    '''
    def __init__(self, residue):
        Residue.__init__(self, residue.bio_res_class)
        self.coords = [atom.coord for atom in self.atoms]

        self.hydrogen = []
        for atom in residue.atoms:
            if atom.element == 'O':
                self.oxygen = atom
            elif atom.element == 'H':
                self.hydrogen.append(atom)
            atom.set_parent(self)

        assert self.oxygen is not None
        assert len(self.hydrogen) == 2

    def oxygen_coordinates(self):
        return self.oxygen.coord

    def H_O_bond_vector(self, hydrogen):
        return hydrogen.coord - self.oxygen_coordinates()


class ActiveSite():
    ''' Methods pertaining to active site information.
    '''
    def __init__(self, index_file, pdb_file):
        self.index_file = index_file
        self.pdb_file = pdb_file

    def read_pdb(self):
        ''' Get lines from pdb file.
        '''
        with open(self.pdb_file, 'r') as f:
            return tuple(line for line in f)

    def read_index(self):
        ''' Get AS line numbers from index file.
        '''
        with open(self.index_file, 'r') as f:
            # there is only 1 line in index files
            index_file_lines = [line for line in f][0]
            return parse.str2list(index_file_lines)

    def get_residues(self, residues):
        ''' Get residues belonging to the active site.
        '''
        line_numbers = self.read_index()
        pdb_lines = self.read_pdb()
        as_lines = [parse.pdbline(pdb_lines[i]) for i in line_numbers]
        # a[3] is residue name, a[5] is residue id
        as_res_in_pdb = {(a[3], int(a[5])) for a in as_lines}
        self.as_res = [r for r in residues if (r.name, r.id) in as_res_in_pdb]

    def get_coords(self):
        ''' Get coordinates of all atoms in the active site.
        '''
        self.coords = []
        for r in self.as_res:
            for a in r.atoms:
                self.coords.append(a.coord)

    def get_central_point(self):
        ''' Get central point of coords of active site.
        '''
        self.get_coords()
        return np.mean(self.as_coords, axis=0)


class ResiduesInFile:
    ''' Container for residues in file. Primary operation is to extract
    all residue data and obtain coordinates and charges of atoms in each
    residue.
    '''
    def __init__(self, pdb, pdb_code, rtp_file):
        self.pdb = pdb
        self.pdb_code = pdb_code
        self.rtp_file = rtp_file

    def extract_residues(self):
        ''' Get all residues in the file.
        '''
        print("Reading file {}".format(self.pdb))
        parser = PDBParser(QUIET=True)
        pdb_file = preprocess_pdb(self.pdb, self.pdb_code, parser)
        struct = parser.get_structure(self.pdb_code, pdb_file)
        return [Residue(r) for r in struct.get_residues()]

    def match_atoms_and_bonds(self, residues, rtp_file):
        ''' Read rtp (main unless specified) and get rtp data for each
        residue in residues.
        '''
        rtp_atoms, rtp_bonds = fileProcessing.RTP.fetch_rtp(rtp_file)
        for residue in residues:
            residue.get_dict_values(rtp_atoms, rtp_bonds, residue.name)
            if not residue.in_dict:
                alt_names = [k for k in rtp_atoms if residue.name in k]
                for name in alt_names:
                    residue.get_dict_values(rtp_atoms, rtp_bonds, name)
                    if residue.in_dict:
                        break
            if not residue.in_dict:
                residue.get_dict_values(rtp_atoms,
                                        rtp_bonds,
                                        residue.name,
                                        forgiving=True)
                if residue.in_dict:
                    w = '''RTP entry for residue {} {} has a different 
                           number of atoms than the model residue. 
                           Accepting this residue anyway.'''
                    warnings.warn(w.format(residue.name, residue.id))
        residues_not_in_dict = [r for r in residues if not r.in_dict]
        for r in residues_not_in_dict:
            w = 'Residue {} {} not found in rtp file.'
            warnings.warn(w.format(r.name, r.id))
            
        return residues

    def process_solvent(self, solvent):
        ''' Differs from process_solute (below) in that histidine
        is not replaced. Other differences may be added.
        '''
        return self.match_atoms_and_bonds(solvent, self.rtp_file)

    def process_solute(self, solute):
        ''' Replace histidine name to match RTP file and get
        atoms and bonds belonging to residues in the RTP.
        '''
        def replace_his(residue):
            ''' HIS is not in RTP file. Keep res, change name.
            '''
            if residue.name == 'HIS':
                residue.name = 'HIE'
            return residue

        for residue in solute:
            residue = replace_his(residue)

        return self.match_atoms_and_bonds(solute, self.rtp_file)

    def separate_components(self, residues):
        ''' Split into solvent (WATER) and solute residues.
        '''
        solvent, solute = [], []
        for residue in residues:
            if residue.is_not_water:
                solute.append(residue)
            else:
                solvent.append(Water(residue))
        return solvent, solute

    def get_residues(self):
        ''' Get all residues in file.
        '''
        all_residues = self.extract_residues()
        solvent, solute = self.separate_components(all_residues)
        solute = self.process_solute(solute)
        self.bad_residues = [r for r in solute if not r.is_good_residue()]
        self.solute = [r for r in solute if r.is_good_residue()]
        self.solvent = self.process_solvent(solvent)
        self.residues = self.solute + self.solvent
        self.solute_coords = [a.coord for r in self.solute for a in r.atoms]
        self.solvent_coords = [a.coord for r in self.solvent for a in r.atoms]

    def get_charges_and_coords(self, residues):
        ''' Extract charge and coords of each atom in the file
        of a certain type (solute/solvent/active site, etc..).
        '''
        unassigned_warning = "Atom{} belonging to " + \
                             "residue {} {} has no assigned charge"

        charges_and_coords = []
        for residue in residues:
            # get H bonded to non-H
            residue.get_bond_partners()
            # rename H in pdb to H in RTP
            residue.sort_hydrogen()
            # get bond partners, now with new H names
            residue.get_bond_partners()
            # get charges from RTP, now that H have proper names
            residue.get_atom_charges()
            # add residue atom info to charges_and_coords
            for atom in residue.atoms:
                cnc = (atom.charge, atom.coord)
                charges_and_coords.append(cnc)

                if atom.charge is None:
                    #print(atom,residue.atom_lines, residue.atoms)
                    continue
                    '''
                    warnings.warn(
                        unassigned_warning.format(
                                atom.fullname,
                                residue.name,
                                residue.id
                            )
                        )'''

        return [c for c in charges_and_coords if c[0] != None]
