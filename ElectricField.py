from spatial.cartesian import *
import Bio

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
    ''' Get the electrostatic field at a given point from the sum
    of other point charges in the model. Result is in units given
    by PDB and is converted after returning for sake of runtime.
    '''
    # separate charges and coords
    charges, coords = zip(*charges_coords)

    # Get cartesian vectors between point charge and target
    vectors = tuple((target - c) for c in iter(coords))
    dists = tuple(magnitude(v) for v in iter(vectors))

    # get the magnitude and i,j,k components of field from point q
    n_c = range(len(coords))
    e_components = ((charges[i] * vectors[i]) / dists[i]**3 for i in n_c)

    # magnitude of electric field on target summed over all contributions,
    return magnitude(sum(e_components))


def get_avg_e_field(targets, charges_coords, cutoff=100):
    ''' Get the average E field felt by a collection of points.
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
            
    field_sum = 0
    for atom in iter(targets):
        ac = atom.coord
        loc = [c for c in charges_coords if dist_le_cutoff(ac, c[1], cutoff)]
        het = tuple(c for c in loc if keep(c[1], atom.exclusions))
        field_sum += electric_field(target=ac, charges_coords=het)
        
    return CONVERSION * field_sum / len(targets)
