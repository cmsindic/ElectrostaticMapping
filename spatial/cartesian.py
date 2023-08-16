import numpy as np
from itertools import product


def magnitude(x):
    ''' Magnitude of an i,j,k vector.
    '''
    return sum(x**2)**0.5


def dist(x, y):
    ''' Euclidian distance between cartesian coordinates.
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
    if abs(z1-z2) > cutoff:
        return False
    if dist(p1, p2) < cutoff:
        return True
    else:
        return False

def target_near_coords(target, coords, cutoff):
    ''' Check if target coord is within cutoff of
    any of coords.
    '''
    for c in coords:
        if dist_le_cutoff(target, c, cutoff):
            return True
    else:
        return False
        
        
class Cube:
    def __init__(self, p1, p2, edge_length):
        self.p1 = np.array(p1)
        self.p2 = np.array(p2)
        self.edge_length = edge_length

    def check_if_contains_point(self, coords):
        print(self.p1,self.p2)
        def contains_point(point):
            for v in iter(point - self.p1):
                if v > self.edge_length or v < 0:
                    return False
            else:
                return True

        for point in iter(coords):
            self.contains_point = contains_point(point)
            if self.contains_point:
                break

    def get_central_point(self):
        return [v + self.edge_length / 2 for v in self.p1]


def generate_grid(coords, res, padding):
    ''' Create a grid of evenly spaced points in the model. Res is
    the number of points per angstrom and padding is the number of
    angstroms to add to each wall of the box.
    '''

    def min_max(arr):
        ''' Return min and max of each array in *args.
        '''
        padded_min = round(min(arr) - padding)
        padded_max = round(max(arr) + padding)
        return [padded_min, padded_max]

    # [(x1, y1, z1), ... (xn, yn, zn)] -->
    # [[x1, ... xn], [y1, ... yn], [z1, ... zn]]
    zipped = tuple(zip(*coords))

    # bounds of box
    min_x, max_x = min_max(zipped[0])
    min_y, max_y = min_max(zipped[1])
    min_z, max_z = min_max(zipped[2])

    # numbers of points in each direction
    range_i = tuple(int(x * res) for x in range((max_x - min_x)))
    range_j = tuple(int(x * res) for x in range((max_y - min_y)))
    range_k = tuple(int(x * res) for x in range((max_z - min_z)))

    # coord values in each axis of grid
    X = [i / res + min_x for i in range_i]
    Y = [i / res + min_y for i in range_j]
    Z = [i / res + min_z for i in range_k]

    # generate grid as dict of dict of dicts
    # with keys as indices so that
    # grid[0][0][0] == min_x, min_y, min_z &
    # grid[-1][-1][-1] == max_x, max_y, max_z
    # Sorry, reader.
    grid = {
            i:{
                j:{
                    k:[X[i], Y[j], Z[k]]
                    for k in range_k
                    }
                for j in range_j
                }
            for i in range_i
        }

    return grid


def generate_cubes(coords, res, padding=4):
    ''' Convert a grid, represented as a dict of dicts
    of dicts of coords, to a series of cube objects
    '''

    def get_range(key):
        ''' Get index/key ranges for dicts with consecutive
        integers as keys. range max = real range max - 1
        for reasons evident in operations below.
        '''
        return range(max(key) - 1)

    grid = generate_grid(coords, res, padding)

    x_keys = grid.keys()
    y_keys = grid[0].keys()
    z_keys = grid[0][0].keys()
    xyz_keys = (x_keys, y_keys, z_keys)
    key_ranges = (get_range(k) for k in xyz_keys)

    cubes = []
    for i, j, k in product(*key_ranges):
        cube = Cube(p1=grid[i][j][k],
                    p2=grid[i+1][j+1][k+1],
                    edge_length=1/res)
        cubes.append(cube)
    return cubes
