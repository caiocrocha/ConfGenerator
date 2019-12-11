import numpy as np
from class_molecule import Molecule

#self.topology.num_bonds = int(line.split()[0])
#    nbonds = 2*self.topology.num_bonds
#    self.topology.bond_list = np.zeros(nbonds, dtype='uint8')
#    count = 0
#    while line:
#        line = inf.readline().split()
#        for i in line:
#            if count < nbonds:
#                self.topology.bond_list[count] = i
#                count += 1

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

def recursive_rotate_list(molecule, rlist, c, start, b):
    if(c < len(rlist) < molecule.num_atoms):
        for i in range(start, molecule.topology.num_bonds*2, 2):
            p = molecule.topology.bond_list[i]
            q = molecule.topology.bond_list[i+1]
            e = rlist[c]
            if(e == p and q != b):
                rlist.append(q)
                start = i+2
                break
            elif(e == q and p != b): 
                rlist.append(p)
                start = i+2
                break
        for i in range(start, molecule.topology.num_bonds*2, 2):
            p = molecule.topology.bond_list[i]
            q = molecule.topology.bond_list[i+1]
            e = rlist[c]
            if(e == p and q != b):
                rlist.append(q)
            elif(e == q and p != b): 
                rlist.append(p)
        recursive_rotate_list(molecule, rlist, c+1, start, b)

def rotate_list(molecule, rlist, a, b):
    for i in range(0, molecule.topology.num_bonds*2, 2):
        p = molecule.topology.bond_list[i]
        q = molecule.topology.bond_list[i+1]
        if(p == b and q != a):
            rlist.append(q)
            start = i+2
            break
        elif(q == b and p != a): 
            rlist.append(p)
            start = i+2
            break
    for i in range(start, molecule.topology.num_bonds*2, 2):
        p = molecule.topology.bond_list[i]
        q = molecule.topology.bond_list[i+1]
        if(p == b and q != a):
            rlist.append(q)
        elif(q == b and p != a): 
            rlist.append(p)
    recursive_rotate_list(molecule, rlist, 0, start, b)
       
def rotate(molecule, a, b, theta, ntimes, filename):   # bond a-b, rotate chain from b // theta in radians
    rlist = []  # list of atoms to be rotated
    rotate_list(molecule, rlist, a, b)
    rlist = list(set(rlist)) # removes duplicates
    
#    For the construction of the rotation matrix around the axis defined by the bond that
#    connects atoms a and b (vector ab, a to b), there needs to be determined ab's unit vector, 
#    which is mathematicaly defined by its original components divided by the square root
#    of its magnitude.
    molecule.write_pdb(filename, 'w+', 0)
    init_rotation_axis(molecule, a, b, theta, ntimes)
    if a < b:
        for i in range(1, ntimes+1):
            for atom in rlist:
                rotate_atom(molecule, atom-1)
            molecule.write_pdb(filename, 'a', i)
    
def calculate_pot_energy(molecule):
    Vb = 0
    for i in range(0, molecule.topology.num_bonds*2, 2):
        a = molecule.topology.bond_list[i]-1
        b = molecule.topology.bond_list[i+1]-1
        btype = ('{}-{}'.format(molecule.topology.type[a], 
                molecule.topology.type[b])) # bond type (e.g. 'c3-c3' or 'c3-hc')
        try:
            K = molecule.topology.bond_types[btype][0]  # elastic constant of the bond
            b0 = molecule.topology.bond_types[btype][1] # equilibrium bond length
        except:
            btype = ('{}-{}'.format(molecule.topology.type[b], 
                    molecule.topology.type[a]))
            try:
                K = molecule.topology.bond_types[btype][0]
                b0 = molecule.topology.bond_types[btype][1]
            except: pass
        b = (((molecule.x[a] - molecule.x[b])**2) + 
             ((molecule.y[a] - molecule.y[b])**2) + 
             ((molecule.z[a] - molecule.z[b])**2))**0.5
        d = b-b0
        Vb += K*(d**2) # acumulates the potential energy of the bond
    return Vb

molecule = Molecule()

#pdb, psf, out_pdb, out_psf = molecule.get_cmd_line()

molecule.read_psf("butane.psf")
molecule.read_pdb("butane_opt.pdb")
molecule.read_frcmod("butane.frcmod")
#molecule.write_pdb("butane_opt_rotated.pdb", 'w+', 0)
#rotate(molecule, 2, 3, 2*np.pi, 360, 'butane_opt_rotated1.pdb')
Vb = calculate_pot_energy(molecule)
