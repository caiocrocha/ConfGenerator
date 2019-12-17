import math
from class_molecule import Molecule

def is_main():
    return __name__ == '__main__'

def init_rotation_axis(molecule, a, b, theta):
    # For the construction of the rotation matrix around the axis defined by the bond that
    # connects atoms a and b (vector ab, a to b), there needs to be determined ab's unit vector, 
    # which is mathematicaly defined by its original components divided by the square root
    # of its magnitude.
    global Ax, Ay, Az, ux, uy, uz, cos0, sin0, t
    Ax = molecule.x[a-1]
    Ay = molecule.y[a-1]
    Az = molecule.z[a-1]
    Bx = molecule.x[b-1]
    By = molecule.y[b-1]
    Bz = molecule.z[b-1]
    ABx = Bx - Ax
    ABy = By - Ay
    ABz = Bz - Az
    M = (ABx*ABx + ABy*ABy + ABz*ABz)**0.5
    # M is the magnitude of the unit vector u
    ux = (Bx - Ax)/M # x component of u
    uy = (By - Ay)/M # y component of u
    uz = (Bz - Az)/M # z component of u
    cos0 = math.cos(theta) # cosine of the angle of rotation
    sin0 = math.sin(theta) # sine of the angle of rotation
    t = 1 - cos0

def rotate_atom(molecule, i):
    x1 = molecule.x[i] - Ax
    y1 = molecule.y[i] - Ay
    z1 = molecule.z[i] - Az
    molecule.x[i] = Ax + (cos0 + t*ux*ux)*x1 + (t*ux*uy - sin0*uz)*y1 + (t*ux*uz + sin0*uy)*z1
    molecule.y[i] = Ay + (t*ux*uy + sin0*uz)*x1 + (cos0 + t*uy*uy)*y1 + (t*uy*uz - sin0*ux)*z1
    molecule.z[i] = Az + (t*ux*uz - sin0*uy)*x1 + (t*uy*uz + sin0*ux)*y1 + (cos0 + t*uz*uz)*z1

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
    init_rotation_axis(molecule, a, b, theta)
    recursive_fast_rotate(molecule, atom=b-1, previous=a-1)
    
def recursive_rotation_list(molecule, atom, previous, rlist):
    k = int(atom*(atom-1)/2)
    for i in range(atom-1, -1, -1):
        if(molecule.topology.bond_matrix[k+i] and i != previous):
            rlist.append(i)
            recursive_rotation_list(molecule, i, atom, rlist)
    for i in range(atom+1, molecule.num_atoms):
        if(molecule.topology.bond_matrix[int(i*(i-1)/2 + atom)] and i != previous):
            rlist.append(i)
            recursive_rotation_list(molecule, i, atom, rlist)
    
def rotate(molecule, rotation_list, a, b, theta, ntimes, write_pdb=False, pdb_name=None):
    init_rotation_axis(molecule, a, b, theta)
    if not write_pdb:
        for i in range(ntimes):
            for atom in rotation_list:
                rotate_atom(molecule, atom)
    else:
        molecule.write_pdb(pdb_name, 'w+', 0)
        for i in range(ntimes):
            for atom in rotation_list:
                rotate_atom(molecule, atom)
            molecule.write_pdb(pdb_name, 'a', i+1)
                      
def get_rotation_list(molecule, a, b):
    rotation_list = []
    recursive_rotation_list(molecule, b-1, a-1, rotation_list)
    return rotation_list

##############################################################################
if is_main():
    
    molecule = Molecule()
    
    # pdb, psf, out_pdb, out_psf = molecule.get_cmd_line()
    
    molecule.read_psf('butane.psf')
    molecule.read_pdb('butane_opt.pdb')
    rotation_list_C3 = get_rotation_list(molecule, 2, 3)
    rotate(molecule, rotation_list_C3, 2, 3, math.pi/180, 360, 
            write_pdb=True, pdb_name='butane_opt_rotated.pdb')
