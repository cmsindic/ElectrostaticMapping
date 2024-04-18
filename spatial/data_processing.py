import numpy as np


def magnitude(x):
    """Calculate the magnitude in 3 dimensions of a len=3 array of distances."""
    return sum(x**2)**0.5

class Project:
    def __init__(self, u):
        """Initialize the Project class with a vector 'u'."""
        self.u = u
    def onto(self, v):
        """Project the vector 'u' onto vector 'v'."""
        v_norm = magnitude(v)
        return (np.dot(self.u, v) / v_norm**2) * v

class AngleBetween:
    def __init__(self, u):
        """Initialize the AngleBetween class with a vector 'u'."""
        self.u = u
    def _and(self, v):
        """Calculate the angle between vectors 'u' and 'v' in degrees."""
        v_norm = magnitude(v)
        u_norm = magnitude(self.u)
        cos_theta = np.dot(self.u, v) / (v_norm * u_norm)
        return np.degrees(np.arccos(cos_theta))


def get_bond_vector(target_coords, other_atom_coords):
    """
    Get the vector of the bond from hydrogen to parent oxygen.
    """
    bond_vector = np.array(target_atom_coords) - np.array(atom_coords)
    return bond_vector


