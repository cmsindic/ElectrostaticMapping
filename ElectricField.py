from spatial.cartesian import *
from scipy.spatial import distance_matrix
import Bio
import numpy as np

# coulomb constant
COUL_CONST = 9 * 10**9
# coulombs per elementary charge
COUL_PER_E = 1.602 * 10**-19
# convert angstroms to meters
ANGSTROM_TO_METER = 10**-10
# convert volts/meter to a.u.
VM_TO_AU = 5.14 * 10**11
# convert E field units from pdb to a.u.
CONVERSION = COUL_CONST * COUL_PER_E / \
             (VM_TO_AU * (ANGSTROM_TO_METER**2))


def electric_field(target, charges_coords):
    ''' Get the electrostatic field felt by a water hydrogen
    atom projected onto bond vector with its parent oxygen.
    '''
    
    # Coordinates of target atom
    target_atom_coords = target.coord
    
    # Bond vector with parent oxygen
    b = target.bond_to_oxygen
    
    # separate charges and coords
    charges, coords = zip(*charges_coords)

    # Get cartesian vectors between point charge and target
    vectors = tuple((target_atom_coords - c) for c in iter(coords))
    dists = tuple(magnitude(v) for v in iter(vectors))

    # get the magnitude and i,j,k components of field from point q
    n_c = range(len(coords))
    e_components = ((charges[i] * vectors[i]) / dists[i]**3 for i in n_c)
    
    net_field = sum(e_components)
    #net_projected_field_mag = MagnitudeOfVector(net_field).projected_onto(b)

    # magnitude of electric field on target summed over all contributions,
    return CONVERSION * magnitude(net_field) #net_projected_field_mag


def get_e_field_list(targets, charges_coords, cutoff):
    ''' Return a list of the fields felt by all water protons 
    in a given partition of water. 
    '''
      
    def keep(coord, exclusions):
        ''' Very explicit way to see if coord is in
        exclusions.
        '''
        x1, y1, z1 = coord
        for point in exclusions:
            x2, y2, z2 = point
            if (x1, y1, z1) == (x2, y2, z2):
                return False
        else:
            return True
            
    fields = []   
    for atom in iter(targets):
        ac = atom.coord
        _, coords = zip(*charges_coords)
        d = enumerate(distance_matrix([ac], coords)[0])
        loc = [charges_coords[i] for i, v in d if v < cutoff]
        het = tuple(c for c in loc if keep(c[1], atom.exclusions))
        fields.append(electric_field(target=atom, charges_coords=het))
    
    return fields


def get_avg_e_field(targets, charges_coords, cutoff=100):
    ''' Get the average E field felt by a collection of points.
    In the main case, this collection is a partition of water 
    protons. 
    '''
    fields = get_e_field_list(targets, charges_coords, cutoff)
    return np.mean(fields)
