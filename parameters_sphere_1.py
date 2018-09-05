import numpy as np

# All numbers in cm
dipole_loc = 7.8
brain_rad = 7.9

sigma_brain = 1. / 300.  # S / cm

# from gmsh sphere_4.geo
brainvol = 32


# measument points
# theta = np.arange(0, 180)
# phi_angle = 0 # -90 to 90

theta, phi_angle = np.mgrid[0:180:1, -90:90:1]
theta = theta.flatten()
phi_angle = phi_angle.flatten()

theta_r = np.deg2rad(theta)
phi_angle_r = np.deg2rad(phi_angle)

rad_tol = 1e-2
x_points = (brain_rad - rad_tol) * np.sin(theta_r) * np.cos(phi_angle_r)
y_points = (brain_rad - rad_tol) * np.sin(theta_r) * np.sin(phi_angle_r)
z_points = (brain_rad - rad_tol) * np.cos(theta_r)

#ele_coords = np.vstack((x_points, y_points, z_points)).T
coords = open('coordinates_NW_later.txt', 'r')
ele_location = coords.readlines()
coords.close()
coordinates = np.zeros([len(ele_location), 3])
for i in range(len(ele_location)):
    crd = ele_location[i].split('\t')
    coordinates[i, :] = crd[0:3]
    for j in range(3):
        coordinates[i, j] = float(coordinates[i, j])
    coordinates[i] = list(coordinates[i])
ele_coords = coordinates[42:75, :]/10

# dipole location - Radial
rad_dipole = {'src_pos': [8., 10., 10.85],
              'snk_pos': [8., 11., 10.75],
              'name': 'rad'}

# # dipole location - Tangential
tan_dipole = {'src_pos': [8., 10.05, 7.8],
              'snk_pos': [8., 13.05, 7.8],
              'name': 'tan'}

# # # dipole location - Mix
mix_dipole = {'src_pos': [8., 8.0353, 7.835],
              'snk_pos': [8., 8.0353, 7.764],
              'name': 'mix'}

dipole_list = [rad_dipole, tan_dipole, mix_dipole]
