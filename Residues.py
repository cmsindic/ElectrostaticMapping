import os
import Bio
import requests
import warnings
import fileProcessing.RTP
import fileProcessing.parse as parse
from fileProcessing.PDB import preprocess_pdb
from spatial.cartesian import *
from textwrap import dedent
from lxml import etree
from io import StringIO
from Bio.PDB import *



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
    
    
    def elements_match_rtp(self, rtp_atoms_entry):
        # get elements of the residue's atoms
        elements = [atom.element for atom in self.atoms].copy()
        # first char in RTP atom name is presumably the element 
        rtp_elements = [a[0][0] for a in rtp_atoms_entry.copy()]
        # dict counting the frequency of each element in the residue
        frequency_dict = {a:elements.count(a) for a in set(elements)}
        # make sure there are at least the same number of unique elements
        assert len(set(rtp_elements) - frequency_dict.keys()) == 0
        # remove a counter for each element as encountered in rtp entry
        for e in rtp_elements:
            frequency_dict[e] -= 1
        for k, v in frequency_dict.items():
            if v != 0: return False
        return True
            
                    
    def search_rtp_dict(self, rtp_atoms_dict, rtp_bonds_dict, name):
        ''' See if residue is in dicts of rtp atom lines and rtp bond lines.
        '''
        assert rtp_atoms_dict.keys() == rtp_bonds_dict.keys()  
        
        if not name in rtp_atoms_dict.keys():
            self.in_dict = False
            return   
            
        else:
            rtp_entry = rtp_atoms_dict[name]
            same_atom_count = len(self.atoms) == len(rtp_entry)
            if same_atom_count and self.elements_match_rtp(rtp_entry):
                self.atom_lines = rtp_atoms_dict[name]
                self.bond_lines = rtp_bonds_dict[name]
                self.in_dict = True
                return     
            else:
                self.in_dict = False
                return

    
    def attempt_dict_match(self, rtp_atoms, rtp_bonds):
        ''' Try to find the residue in the dicts of atoms and bonds
        created from lines the the rtp file. Identify if PDB residue
        is in both dictionaries. 
        '''
        
        # 'HIS' is generally not in RTP file, try 'HIE'
        #if self.name == 'HIS': self.name = 'HIE'
        
        # see if entirety of residue data is found in rtp file
        self.search_rtp_dict(rtp_atoms, rtp_bonds, self.name)
        
        # if not, the residue may be a variation
        if not self.in_dict:
            
            # variations of HIS are specially named in the RTP file
            if self.name == 'HIS':
                alt_names = []
                for n in ('HIE','HID','HIP'):
                    alt_names += [x for x in rtp_atoms.keys() if n in x]
            else:
                # e.g. ASN may be CASN in the rtp file ('ASN' is in 'CASN')
                alt_names = (k for k in rtp_atoms.keys() if self.name in k)
                
            # search the RTP dicts for any of the new alt_names
            for name in alt_names:
                self.search_rtp_dict(rtp_atoms, rtp_bonds, name)
                if self.in_dict: break


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
            # add H bonded to atom to H belonging to res
            self.rtp_hydrogen += atom.rtp_hydrogen
            # n H bonded to atom
            atom.n_children = len(atom.rtp_hydrogen) 
            # init empty array of H assigned to atom
            atom.bonded_h = [] 

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
            for i, h in enumerate(atom.bonded_h):
                h.name = atom.rtp_hydrogen[i]


    def get_atom_charges(self):
        ''' Get atomic partial charge for a given atom.
        '''
        
        def print_guess_message(resname, res_id,
                                atom_name, rtp_name):
            ''' Warning to notify user of atoms that are guessed.
            '''
            m = dedent("For residue {} {}: " + "guessing that " + 
                       "model atom {} is atom {} " + "in the RTP file")      
            print(m.format(resname, res_id, atom_name, rtp_name))
        
        def assign_charge(atom, op, rtp_atoms):
            ''' Assign atom charge, rename atom, and
            and update rtp_atoms by remove the atom upon
            matching charges/names.
            '''                              
            # exact match between rtp name (x[0]) and model atom name
            if op == 'EQ':
                match = lambda x: x[0] == atom.name  
            # match between 1st char of rtp name & 1st char of model atom name
            elif op == 'FIRST_CHAR':
                match = lambda x: x[0][0] == atom.name[0]
                
            # recall: line = name_in_rtp, _, charge, _
            matches = [line for line in rtp_atoms if match(line)]
            
            # if there are matches between model atom and rtp atom
            if len(matches) > 0:
                # select the first 
                rtp_name, _, charge, _ = matches[0] 
                # print notification of atom guess if not exact match
                if op == 'FIRST_CHAR':
                    print_guess_message(self.name, self.id,
                                        atom.name, rtp_name)
                # update atom name and charge
                atom.name, atom.charge = rtp_name, charge 
                # remove the matching rtp line from candidates
                rtp_atoms.remove(matches[0])
            return atom, rtp_atoms

        # assign_charge will gradually remove matched atoms
        # from self.atom_lines as they are assigned, so give
        # a copy as an argument instead
        rtp_atoms = self.atom_lines.copy()
        
        # try to find atoms in RTP that exactly match names in model
        for atom in self.atoms:
            atom.charge = None
            atom, rtp_atoms = assign_charge(atom, 'EQ', rtp_atoms)
        
        # this will rename atoms with partially matching names
        chargeless_atoms = [a for a in self.atoms if a.charge==None]
        for atom in chargeless_atoms:
            atom, rtp_atoms = assign_charge(atom, 'FIRST_CHAR', rtp_atoms)
    
    
    def check_for_chargeless_atoms(self):
        ''' Warn user about atoms that are have not been assigned charges.
        '''
        w = "Atom {} belonging to residue {} {} has no assigned charge"
        chargeless_atoms = [a for a in self.atoms if a.charge == None]
        for atom in chargeless_atoms:        
            warnings.warn(w.format(atom.fullname, self.name, self.id))
          
          
    def charges(self):
        ''' Yield charges of atoms in residue.
        '''
        for atom in self.atoms:
            yield atom.charge
         
         
    def coords(self):
        ''' Yield coordinates of atoms in residue.
        '''
        for atom in self.atoms:
            yield atom.coord
    
    
    def get_water_attributes(self):
        ''' Get properties unique to water, for the identification 
        of oxygen and hydrogen in each water molecule.
        '''
        for atom in self.atoms:
            self.hydrogen = []
            if atom.element == 'O':
                self.oxygen = atom
            elif 'H' in atom.name:
                self.hydrogen.append(atom)
            atom.exclusions = tuple(a.coord for a in self.atoms)

"""
class Water(Residue):
    ''' Special methods for handling water if electrostatic fields
    are to be measured as they are felt by water's components.
    '''
    def __init__(self, residue):
        Residue.__init__(self, residue.bio_res_class)
        self.hydrogen = []
        for atom in residue.atoms:
            if atom.element == 'O':
                self.oxygen = atom
                self.oxygen.coord = atom.coord
            elif atom.element == 'H':
                self.hydrogen.append(atom)
            atom.set_parent(self)
        
        # Water must have 1 oxygen and 2 hydrogen
        assert self.oxygen is not None
        assert len(self.hydrogen) == 2


    def H_O_bond_vector(self, hydrogen):
        ''' Get the bond vector between a hydrogen and its 
        parent oxygen atom. Important for projecting electric
        field along this bond vector. See ElectricField.py.
        '''
        return hydrogen.coord - self.oxygen.coord
"""


class ActiveSite():
    ''' Methods pertaining to active site information.
    '''
    def __init__(self, pdb_id):
        self.pdb_id = pdb_id


    def get_rcsb_tree(self):
        ''' Return the xml tree of the rscb page.
        '''
        parser = etree.HTMLParser()
        base_rcsb_url = 'https://www.rcsb.org/structure/'
        url = base_rcsb_url + self.pdb_id
        page = requests.get(url)
        html = page.content.decode("utf-8")
        return etree.parse(StringIO(html), parser=parser)


    def fetch_uniprot_id(self):
        ''' Resolve the identity of the uniprot object linked 
        to the pdb model by extracting the name from the uniprot 
        url that is linked in the rscb page for the model.
        '''
        tree = self.get_rcsb_tree()
        # Get anchor tags
        refs = tree.xpath("//a")
        # Get refs
        links = [link.get('href', '') for link in refs]
        # Isolate uniprot url links
        uniprot_links = [l for l in links if 'uniprot.org' in l]
        # Assert there is only 1 uniprot link
        assert len(uniprot_links)==1
        # Ex: 'https://www.uniprot.org/uniprot/P00942' --> 'P00942'
        return uniprot_links[0].split('/')[-1]
        
        
    def fetch_uniprot_text(self):
        ''' Get the raw text of the uniprot text file corresponding
        with the pdb model.
        '''
        uniprot_base_url = 'https://rest.uniprot.org/uniprotkb/'
        uniprot_id = self.fetch_uniprot_id()
        url = uniprot_base_url + uniprot_id + '.txt'
        page = requests.get(url)
        return page.text
        
        
    def resolve_res_ids(self):
        ''' Extract the residue ID numbers of the active site
        residues from the uniprot text file for the pdb model.
        '''
        # Text file from uniprot detailing model features
        raw_text = self.fetch_uniprot_text()
        # Split text into list of lines
        text = raw_text.split('\n')
        # Split each line into list of words
        parsed = (parse.text_to_list(line) for line in text)
        # Isolate lines containing active site info
        act_site_lines = (line for line in parsed if 'ACT_SITE' in line)
        # Residue IDs are the last words in the lines
        return tuple(int(line[-1]) for line in act_site_lines)
        
                
    def get_residues(self, residues):
        ''' Get residues belonging to the active site.
        '''
        act_site_ids = self.resolve_res_ids()
        self.as_res = [r for r in residues if r.id in act_site_ids]

        print("Processing active site")
        print("Active site residues: ")
        
        #### Ask if assumed active site residues are true ####
        
        # Because thare are often multiple chains, make a 
        # unique list of residue names/ids to inquire about
        no_dup = list(set((r.name, r.id) for r in self.as_res))
        
        # Ask if each residue is valid
        for res in no_dup:
            print(res[0], res[1])
            keep = bool(input("Keep this residue? (y/n) \n"))
            if not keep:
                no_dup.remove(res)
        
        # Residues will now be removed from multiple chains,
        # even though user is asked one the basis of name & id
        self.as_res = [r for r in self.as_res if (r.name, r.id) in no_dup]


    def get_coords(self):
        ''' Get coordinates of all atoms in the active site.
        '''
        self.coords = [a.coord for r in self.as_res for a in r.atoms]
        
        
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
    def __init__(self, pdb, pdb_id, rtp_file):
        self.pdb = pdb
        self.pdb_id = pdb_id
        self.rtp_file = rtp_file

    def extract_residues(self):
        ''' Get all residues in the file.
        '''
        print("Reading file {}".format(self.pdb))
        parser = PDBParser(QUIET=True)
        pdb_file = preprocess_pdb(self.pdb, self.pdb_id, parser)
        struct = parser.get_structure(self.pdb_id, pdb_file)
        return [Residue(r) for r in struct.get_residues()]

    def fetch_rtp_data(self, residues, rtp_file):
        ''' Read rtp (main unless specified) and get rtp data for each
        residue in 'residues'.
        '''
        # lines from rtp file containing atom and bond information 
        # for each residue, in dict format with renames as keys
        rtp_atoms, rtp_bonds = fileProcessing.RTP.fetch_rtp(rtp_file)
        for residue in residues:
            residue.attempt_dict_match(rtp_atoms, rtp_bonds)
            
        residues_not_in_dict = [r for r in residues if not r.in_dict]
        for r in residues_not_in_dict:
            w = 'Residue {} {} not found in rtp file.'
            warnings.warn(w.format(r.name, r.id))
        
        return residues

    def process_solvent(self, solvent):
        ''' Get atoms and bonds belonging to water in the RTP.
        '''
        return self.fetch_rtp_data(solvent, self.rtp_file)

    def process_solute(self, solute):
        ''' Get atoms and bonds belonging to residues in the RTP.
        '''
        return self.fetch_rtp_data(solute, self.rtp_file)

    def partitions(self):
        ''' Split into solvent (WATER) and solute residues.
        '''
        residues = self.extract_residues()
        solute = [r for r in residues if r.is_not_water]
        solvent = [r for r in residues if not r.is_not_water]
        return solvent, solute

    def get_residues(self):
        ''' Get all residues in file.
        '''
        solvent, solute = self.partitions()
        self.solute = [r for r in self.process_solute(solute) if r.in_dict]
        self.solvent = [r for r in self.process_solvent(solvent) if r.in_dict]
        self.solute_coords = [a.coord for r in self.solute for a in r.atoms]
    
    def get_all_molecules(self):
        ''' Return all residues in model.
        '''
        return tuple(self.solute + self.solvent)

    def convert_solvent_to_water_objects(self):
        ''' Add features to solvent unique to water, namely the identifying 
        properties of the parent oxygen coordinate and the bond vectors
        between the H and O of the molecules. 
        '''
        for wat in self.solvent:
            wat.get_water_attributes()
            
    def get_charges_and_coords(self, residues):
        ''' Extract charge and coords of each atom in the file
        of a certain type (solute/solvent/active site, etc..).
        '''
        charges, coords = [], []
        for residue in residues:
            # get H bonded to non-H
            residue.get_bond_partners()
            # rename H in pdb to H in RTP
            residue.sort_hydrogen()
            # get bond partners, now with new H names
            residue.get_bond_partners()
            # get charges from RTP, now that H have proper names
            residue.get_atom_charges()
            # check for atoms that are missing charges
            residue.check_for_chargeless_atoms()
            # add residue atom info to charges and coords
            charges += residue.charges()
            coords += residue.coords()
            
        return tuple(c for c in zip(charges, coords) if c[0] != None)
