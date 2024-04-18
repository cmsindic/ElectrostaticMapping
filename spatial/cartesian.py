import numpy as np
from itertools import product
from scipy.spatial import distance_matrix


def magnitude(x):
    '''Magnitude of an i, j, k vector.
    '''
    return sum(x**2)**0.5


def dist(x, y):
    '''Euclidean distance between Cartesian coordinates.
    '''
    d = 0
    for i in range(3):
        d += (x[i] - y[i])**2
    return d**0.5


def dist_le_cutoff(p1, p2, cutoff):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    if abs(x1 - x2) > cutoff:
        return False
    if abs(y1 - y2) > cutoff:
        return False
    if abs(z1 - z2) > cutoff:
        return False
    if dist(p1, p2) < cutoff:
        return True
    else:
        return False


def target_near_coords(target, coords, cutoff):
    '''Check if the target coordinate is within 
    the cutoff of any of the coordinates.
    '''
    d = distance_matrix([target], coords)[0]
    return any(d < cutoff)


def magnitude(x):
    '''Calculate the magnitude in 3 dimensions of
    a len=3 array of distances.
    '''
    return sum(x**2)**0.5


class AngleBetween:
    '''Calculate the angle between two vectors.
    '''
    def __init__(self, u):
        '''Initialize the AngleBetween class with
        a vector 'u'.
        '''
        self.u = u

    def _and(self, v):
        '''Calculate the angle between vectors 'u'
        and 'v' in degrees.
        '''
        v_norm = magnitude(v)
        u_norm = magnitude(self.u)
        cos_theta = np.dot(self.u, v) / (v_norm * u_norm)
        return np.degrees(np.arccos(cos_theta))


class Project:
    '''Project a vector onto another.
    '''
    def __init__(self, u):
        '''Initialize the Project class with
        a vector 'u'.
        '''
        self.u = u

    def onto(self, v):
        '''Project the vector 'u' onto vector 'v'.
        '''
        v_norm = magnitude(v)
        return (np.dot(self.u, v) / v_norm**2) * v


def get_bond_vector(target_atom_coords, other_atom_coords):
    '''Get the vector of the bond from 
    hydrogen to parent oxygen.
    '''
    a = np.array(target_atom_coords)
    b = np.array(other_atom_coords) 
    return a - b 


class MagnitudeOfVector:
    '''Return the magnitude of the net field 
    vector on a hydrogen projected onto the 
    bond vector with its parent oxygen.
    '''
    def __init__(self, u):
        self.u = u

    def projected_onto(self, v):
        u = self.u
        projected_field = Project(u).onto(v)
        angle = AngleBetween(u)._and(v)
        if 0 < angle <= 90:
            return magnitude(projected_field)
        elif 90 < angle <= 180:
            return -magnitude(projected_field)
