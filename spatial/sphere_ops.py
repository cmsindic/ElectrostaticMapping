import numpy as np
import CoordSphere as CoordSphere

def create_sphere(center_point,radius, resolution):
    ''' Create linspaces in x,y,z representing a shpere around
    center_point with radius and resolution being number of samples
    in the linspace.
    '''

    def convert_linspace_to_cartesian(x, y, z):
        ''' Convert 3 linespaces -- x,y,z -- to cartesian coordinates.
        '''
        coords = []
        for i in range(len(x)):
            for j in range(len(x[i])):
                coords.append((x[i][j], y[i][j], z[i][j]))
        return coords

    cx, cy, cz = center_point
    phi = np.linspace(0, 2*np.pi, 2*resolution)
    theta = np.linspace(0, np.pi, resolution)
    theta, phi = np.meshgrid(theta, phi)
    r_xy = radius*np.sin(theta)
    x = cx + np.cos(phi) * r_xy
    y = cy + np.sin(phi) * r_xy
    z = cz + radius * np.cos(theta)
    return convert_linspace_to_cartesian(x, y, z)

def generate_concentric_spheres(
    radii,
    central_point,
    charge_and_coords,
    consider_physics=True,
    cutoff=2
    ):
    ''' Generate concentric spheres of increasing radii around a
    central point.
    '''
    spheres = {}
    for rad in radii:
        sphere = CoordSphere(central_point, rad)
        if consider_physics:
            sphere.remove_VDW_violations(charge_and_coords, cutoff=cutoff)
        spheres[rad] = sphere.get_average_field_around_sphere(charge_and_coords)

class CoordSphere:
    def __init__(
            self,
            center_point,
            radius,
            resolution=7
        ):
        ''' Create a sphere around center_point with a given radius;
        first as linespace, then convert to cartesian.
        '''
        self.sphere_coords = create_sphere(
            center_point,
            radius,
            resolution)

    def remove_VDW_violations(self, charge_and_coords, cutoff):
        ''' Remove points in sphere coords too close to atoms.
        '''
        sphere_coords_copy = self.sphere_coords
        for s in self.sphere_coords:
            for _, coords in charge_and_coords:
                if dist(s, coords) < cutoff:
                    sphere_coords_copy.remove(s)
                    break
            else:
                continue
        self.sphere_coords = sphere_coords_copy

    def get_average_field_around_sphere(self, charge_and_coords):
        ''' Get the average field felt by points around the sphere.
        '''
        return ElectricField.get_avg_e_field(
            target_coords=self.sphere_coords,
            charge_and_coords=charge_and_coords
        )
