
import os
import numpy as np
from dolfin import *
import parameters_sphere_1 as params

def extract_pots(phi, positions):
    compt_values = np.zeros(positions.shape[0])
    for ii in range(positions.shape[0]):
        compt_values[ii] = phi(positions[ii, :])
    return compt_values

def main_1shell_fem(mesh, subdomains, boundaries, skull_cond, src_pos, snk_pos):
    sigma_B = Constant(params.sigma_brain)
    sigma_A = Constant(0.)

    V = FunctionSpace(mesh, "CG", 2)
    v = TestFunction(V)
    u = TrialFunction(V)

    phi = Function(V)
    dx = Measure("dx")(subdomain_data=subdomains)
    ds = Measure("ds")(subdomain_data=boundaries)
    a = inner(sigma_B * grad(u), grad(v))*dx(params.brainvol)
    

    x_pos, y_pos, z_pos = src_pos
    R = 1.2
    f_src = Expression("exp(-(pow(x[0]- x_pos, 2) + pow(x[1] - y_pos, 2) + pow(x[2] + z_pos, 2)) /2.*R)",
                       x_pos=x_pos, y_pos=y_pos, z_pos=z_pos, R=R, degree=3)
    f_snk = Expression("-exp(-(pow(x[0]- x_pos, 2) + pow(x[1] - y_pos, 2) + pow(x[2] + z_pos, 2)) /2.*R)",
                       x_pos=x_pos, y_pos=y_pos, z_pos=z_pos, R=R, degree=3)
    f = f_src + f_snk
    L = f*v*dx(32)


#    L = Constant(0)*v*dx
    A = assemble(a)
    b = assemble(L)

#    x_pos, y_pos, z_pos = src_pos
#    point = Point(x_pos, y_pos, z_pos)
#    delta = PointSource(V, point, 1.)
#    delta.apply(b)
#
#    x_pos, y_pos, z_pos = snk_pos
#    point1 = Point(x_pos, y_pos, z_pos)
#    delta1 = PointSource(V, point1, -1.)
#    delta1.apply(b)

    bc = DirichletBC(V, Constant(0), boundaries, 30) # Surface of the grnd ele
    bc.apply(A, b)
    
    solver = KrylovSolver("cg", "petsc_amg")
    solver.parameters["maximum_iterations"] = 1000
    solver.parameters["absolute_tolerance"] = 1E-8
    solver.parameters["monitor_convergence"] = True
    solver.parameters["error_on_nonconvergence"] = True
    info(solver.parameters, verbose=True)
    set_log_level(PROGRESS)

    solver.solve(A, phi.vector(), b)

    ele_pos_list = params.ele_coords
    vals = extract_pots(phi, ele_pos_list)
    vtkfile = File('sphere_1_gaussian/gauss_solution.pvd')
    vtkfile << (phi)
#    vtkfile = File('sphere_1_gaussian/pots_values.pvd')
#    vtkfile << (vals)
    np.save(os.path.join('results', 'values.npy'), vals)
    return vals


if __name__ == '__main__':
    print('Loading meshes')
    mesh = Mesh("sphere_1.xml")
    subdomains = MeshFunction('size_t', mesh, "sphere_1_physical_region.xml")
    boundaries = MeshFunction('size_t', mesh, "sphere_1_facet_region.xml")
    for dipole in params.dipole_list:
        print('Now computing FEM for dipole: ', dipole['name'])
        src_pos = dipole['src_pos']
        snk_pos = dipole['snk_pos']
        print('Done loading meshes')
        fem_20 = main_1shell_fem(mesh, subdomains, boundaries,
                                 params.sigma_brain, src_pos, snk_pos)
        print('Done 4Shell-FEM-20')
        f = open(os.path.join('results',
                              'sphere_1_Numerical_' + dipole['name'] + '.npz'), 'w')
        np.savez(f, fem_20=fem_20)
        f.close()
