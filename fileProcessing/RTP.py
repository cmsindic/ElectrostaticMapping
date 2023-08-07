class RTPline():
    ''' To contain information about lines in RTP file'''
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
        self.content = [x for x in line if x not in ('', '[', ']')]
        self.is_empty = len(self.content) == 0

def fetch_rtp(rtp_file):
    '''Get lines with charges, atom names, etc.,
    from rtp_atoms.rtp. Return dicts with resnames
    as keys
    '''
    with open(rtp_file, 'r') as f:
        rtp_file = [line for line in f]

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

            # format line if it is an atom line
            if line_type == 'ATOM':
                id, atom, charge, indx = line.content
                formatted_line = (id, atom, float(charge), int(indx))
                atom_lines[resname].append(formatted_line)
            # return raw line if bond line
            elif line_type == 'BOND':
                bond_lines[resname].append(line.content)

        # Skip all lines added to dicts to look for new res start line
        i += j

    return atom_lines, bond_lines
