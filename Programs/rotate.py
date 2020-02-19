import math

def set_rotation_angle(theta):
    global cos0, sin0, t
    cos0 = math.cos(theta)  # cosine of the angle of rotation
    sin0 = math.sin(theta)  # sine of the angle of rotation
    t = 1 - cos0

def init_rotation_axis(molecule, a, b, theta):
    # For the construction of the rotation matrix around the axis defined by the bond that
    # connects atoms a and b (vector ab, a to b), there needs to be determined ab's unit vector, 
    # which is mathematicaly defined by its original components divided by the square root
    # of its magnitude.
    global Ax, Ay, Az, ux, uy, uz
    Ax = molecule.x[a]
    Ay = molecule.y[a]
    Az = molecule.z[a]
    Bx = molecule.x[b]
    By = molecule.y[b]
    Bz = molecule.z[b]
    ABx = Bx - Ax
    ABy = By - Ay
    ABz = Bz - Az
    M = (ABx*ABx + ABy*ABy + ABz*ABz)**0.5
    # M is the magnitude of the unit vector u
    if M != 0:
        ux = ABx/M # x component of u
        uy = ABy/M # y component of u
        uz = ABz/M # z component of u
    else:
        ux = 0
        uy = 0
        uz = 0
    set_rotation_angle(theta)

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
    init_rotation_axis(molecule, a-1, b-1, theta)
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
    init_rotation_axis(molecule, a-1, b-1, theta)
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
    recursive_rotation_list(molecule, b, a, rotation_list)
    return rotation_list
