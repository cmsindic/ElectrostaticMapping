import numpy as np
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

def electric_field(target, charges_coords, exclude=[]):
    ''' Get the electrostatic field at a given point from the sum
    of other point charges in the model.
    '''
    # if target is an atom, other atoms in molecule must be excluded.
    # remove charges/coords belonging to the same molecule as target
    exclude = [list(x) for x in exclude]
    def exclude_coords(coord): return list(coord[1]) in exclude
    charges_coords = tuple(c for c in charges_coords if not exclude_coords(c))

    # separate charges and coords
    charges, coords = zip(*charges_coords)

    # convert charges from elementary to Coulombic
    charges = tuple(COUL_CONST * COUL_PER_E * c for c in charges)

    # Get cartesian vectors between point charge and target
    vectors = tuple((target - c) * ANGSTROM_TO_METER for c in coords)
    dists = tuple(magnitude(v) for v in vectors)

    # get the magnitude and i,j,k components of field from point q
    n_c = range(len(coords))
    e_mags = tuple(charges[i] / dists[i]**2 for i in n_c)
    e_components = ((e_mags[i] * vectors[i]) / dists[i] for i in n_c)

    # magnitude of electric field on target summed over all contributions,
    # in atomic units
    return sum(e_components) / VM_TO_AU

def get_avg_e_field(targets, charges_coords, exclude=[], cutoff=np.Inf):
    ''' Get the average E field felt by a collection of points.
    '''
    field_sum = 0
    for s in targets:

        if isinstance(s, Bio.PDB.Atom.Atom):
            exclude = [a.coord for a in s.parent.atoms]
            s = s.coord

        if cutoff != np.Inf:
            local_coords = tuple(c for c in charges_coords
                                 if dist(c[1], s) < cutoff)
        else:
            local_coords = charges_coords

        field_sum += electric_field(target=s,
                                    charges_coords=local_coords,
                                    exclude=exclude)

    return field_sum / len(targets)
