import numpy as np
from class_molecule import Molecule

def init_rotation_axis(molecule, a, b, theta, times):
#    For the construction of the rotation matrix around the axis defined by the bond that
#    connects atoms a and b (vector ab, a to b), there needs to be determined ab's unit vector, 
#    which is mathematicaly defined by its original components divided by the square root
#    of its magnitude.
    global Ax, Ay, Az, ux, uy, uz, delta_theta, cos0, sin0, t
    Ax = molecule.x[a-1]
    Ay = molecule.y[a-1]
    Az = molecule.z[a-1]
    Bx = molecule.x[b-1]
    By = molecule.y[b-1]
    Bz = molecule.z[b-1]
    M = (((Bx - Ax)**2) + ((By - Ay)**2) + ((Bz - Az)**2))**0.5
    # M is the magnitude of the unit vector u
    ux = (Bx - Ax)/M # x component of u
    uy = (By - Ay)/M # y component of u
    uz = (Bz - Az)/M # z component of u
    delta_theta = theta/times
    cos0 = np.cos(delta_theta) # cosine of the angle of rotation
    sin0 = np.sin(delta_theta) # sine of the angle of rotation
    t = 1 - cos0
    return delta_theta

def rotate_atom(molecule, i):
    x1 = molecule.x[i] - Ax
    y1 = molecule.y[i] - Ay
    z1 = molecule.z[i] - Az
    molecule.x[i] = Ax + (cos0 + t*(ux**2))*x1 + (t*ux*uy - sin0*uz)*y1 + (t*ux*uz + sin0*uy)*z1
    molecule.y[i] = Ay + (t*ux*uy + sin0*uz)*x1 + (cos0 + t*(uy**2))*y1 + (t*uy*uz - sin0*ux)*z1
    molecule.z[i] = Az + (t*ux*uz - sin0*uy)*x1 + (t*uy*uz + sin0*ux)*y1 + (cos0 + t*(uz**2))*z1

def recursive_fast_rotate(molecule, atom, previous):
    k = int(atom*(atom-1)/2)
    for i in range(atom-1, -1, -1):
        if(molecule.topology.bond_matrix[k+i] and i != previous):
            rotate_atom(molecule, i)
            recursive_fast_rotate(molecule, i, atom)
    for i in range(atom+1, molecule.num_atoms):
        if(molecule.topology.bond_matrix[int(i*(i-1)/2 + atom)] and i != previous):
            rotate_atom(molecule, i)
            recursive_fast_rotate(molecule, i, atom)
            
def fast_rotate(molecule, a, b, theta): # bond a-b, rotate chain from b // theta in radians
    init_rotation_axis(molecule, a, b, theta, times=1)
    recursive_fast_rotate(molecule, atom=b-1, previous=a-1)
    
def recursive_rotate_list(molecule, atom, previous, rlist):
    k = int(atom*(atom-1)/2)
    for i in range(atom-1, -1, -1):
        if(molecule.topology.bond_matrix[k+i] and i != previous):
            rlist.append(i)
            recursive_rotate_list(molecule, i, atom, rlist)
    for i in range(atom+1, molecule.num_atoms):
        if(molecule.topology.bond_matrix[int(i*(i-1)/2 + atom)] and i != previous):
            rlist.append(i)
            recursive_rotate_list(molecule, i, atom, rlist)
    
def rotate(molecule, a, b, theta, ntimes, filename): 
    init_rotation_axis(molecule, a, b, theta, ntimes)
    rlist = []
    molecule.write_pdb(filename, 'w+', 0)
    recursive_rotate_list(molecule, b-1, a-1, rlist)
    for i in range(ntimes):
        for atom in rlist:
            rotate_atom(molecule, atom)
        molecule.write_pdb(filename, 'a', i+1)
                      
molecule = Molecule()

#pdb, psf, out_pdb, out_psf = molecule.get_cmd_line()

molecule.read_psf("butane.psf")
molecule.read_pdb("butane_opt.pdb")
rotate(molecule, 2, 3, 2*np.pi, 360, 'butane_opt.pdb')
#fast_rotate(molecule, a=2, b=3, theta=np.pi/2)
#molecule.write_pdb("butane_opt1.pdb", 'w+', 0)
#theta = np.pi/180
#for i in range(360):
#    fast_rotate(molecule, 2, 3, theta)
#    fast_rotate(molecule, 3, 2, theta)
#    molecule.write_pdb("butane_opt_rotated.pdb", 'a', i+1)
