import numpy as np
from class_molecule import Molecule
import matplotlib.pyplot as plt

#def recursive_rotate_list(molecule, r, c, start, b):
#    if(c < len(r) < molecule.num_atoms):
#        for i in range(start, molecule.topology.num_bonds*2, 2):
#            p = molecule.topology.bond_list[i]
#            q = molecule.topology.bond_list[i+1]
#            e = r[c]
#            if(e == p and q != b):
#                r.append(q)
#                start = i+2
#                break
#            elif(e == q and p != b): 
#                r.append(p)
#                start = i+2
#                break
#        for i in range(start, molecule.topology.num_bonds*2, 2):
#            p = molecule.topology.bond_list[i]
#            q = molecule.topology.bond_list[i+1]
#            e = r[c]
#            if(e == p and q != b):
#                r.append(q)
#            elif(e == q and p != b): 
#                r.append(p)
#        recursive_rotate_list(molecule, r, c+1, start, b)
#
#def rotate_list(molecule, r, a, b):
#    for i in range(0, molecule.topology.num_bonds*2, 2):
#        p = molecule.topology.bond_list[i]
#        q = molecule.topology.bond_list[i+1]
#        if(p == b and q != a):
#            r.append(q)
#            start = i+2
#            break
#        elif(q == b and p != a): 
#            r.append(p)
#            start = i+2
#            break
#    for i in range(start, molecule.topology.num_bonds*2, 2):
#        p = molecule.topology.bond_list[i]
#        q = molecule.topology.bond_list[i+1]
#        if(p == b and q != a):
#            r.append(q)
#        elif(q == b and p != a): 
#            r.append(p)
#    recursive_rotate_list(molecule, r, 0, start, b)
#       
#def rotate(molecule, a, b, theta):   # bond a-b, rotate chain from b // theta in radians
#    r = list()  # list of atoms to be rotated
#    rotate_list(molecule, r, a, b)
#    r = list(set(r)) # removes duplicates
#    
##    For the construction of the rotation matrix around the axis defined by the bond that
##    connects atoms a and b (vector ab, a to b), there needs to be determined ab's unit vector, 
##    which is mathematicaly defined by its original components divided by the square root
##    of its magnitude.
#    
#    Ax = molecule.x[a-1]
#    Ay = molecule.y[a-1]
#    Az = molecule.z[a-1]
#    Bx = molecule.x[b-1]
#    By = molecule.y[b-1]
#    Bz = molecule.z[b-1]
#    M = (((Bx - Ax)**2) + ((By - Ay)**2) + ((Bz - Az)**2))**0.5
#    # M is the magnitude of the unit vector u
#    ux = (Bx - Ax)/M # x component of u
#    uy = (By - Ay)/M # y component of u
#    uz = (Bz - Az)/M # z component of u
#    cos0 = np.cos(theta) # cosine of the angle of rotation
#    sin0 = np.sin(theta) # sine of the angle of rotation
#    t = 1 - cos0
#    for i in range(len(r)):
#        j = r[i]-1
#        x1 = molecule.x[j] - Ax
#        y1 = molecule.y[j] - Ay
#        z1 = molecule.z[j] - Az
#        molecule.x[j] = round(((cos0 + t*(ux**2))*x1 + (t*ux*uy - sin0*uz)*y1 + (t*ux*uz + sin0*uy)*z1), 3) + Ax
#        molecule.y[j] = round(((t*ux*uy + sin0*uz)*x1 + (cos0 + t*(uy**2))*y1 + (t*uy*uz - sin0*ux)*z1), 3) + Ay
#        molecule.z[j] = round(((t*ux*uz - sin0*uy)*x1 + (t*uy*uz + sin0*ux)*y1 + (cos0 + t*(uz**2))*z1), 3) + Az
#        # rotation matrix applied to point
#    
#def calculate_pot_energy(molecule):
#    for i in range(molecule.num_bonds*2):            
#        for j in range(len(molecule.topology.bond_types)):
#            bond = (molecule.topology.bond_types[j].split('-')) # bond type (e.g. 'c3-c3' or 'c3-hc')
#            ind_a = molecule.topology.bond_list[i]-1 # index of atom a
#            ind_b = molecule.topology.bond_list[i+1]-1 # index of atom b
#            atom_a = molecule.topology.type[ind_a] # atom type (e.g. 'c3' or 'hc')
#            atom_b = molecule.topology.type[ind_b]
#            if((atom_a == bond[0] and atom_b == bond[1]) or (atom_a == bond[1] and atom_b == bond[0])):
#                K = molecule.topology.k_elastic[j] # elastic constant of the bond
#                # b is the distance between atoms a and b
#                b0 = molecule.topology.bond_lengths[j] # b0 is the equilibrium bond length
#                break

def init_rotation_axis(molecule, a, b):
#    For the construction of the rotation matrix around the axis defined by the bond that
#    connects atoms a and b (vector ab, a to b), there needs to be determined ab's unit vector, 
#    which is mathematicaly defined by its original components divided by the square root
#    of its magnitude.
    global Ax, Ay, Az, ux, uy, uz
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

def rotate_atom(molecule, i, cos0, sin0, t):
    x1 = molecule.x[i] - Ax
    y1 = molecule.y[i] - Ay
    z1 = molecule.z[i] - Az
    molecule.x[i] = round(((cos0 + t*(ux**2))*x1 + (t*ux*uy - sin0*uz)*y1 + (t*ux*uz + sin0*uy)*z1), 3) + Ax
    molecule.y[i] = round(((t*ux*uy + sin0*uz)*x1 + (cos0 + t*(uy**2))*y1 + (t*uy*uz - sin0*ux)*z1), 3) + Ay
    molecule.z[i] = round(((t*ux*uz - sin0*uy)*x1 + (t*uy*uz + sin0*ux)*y1 + (cos0 + t*(uz**2))*z1), 3) + Az
    
def fast_rotate(molecule, a, b, theta): # bond a-b, rotate chain from b // theta in radians

    def recursive_fast_rotate1(molecule, atom):
        for i in range(atom+1, molecule.num_atoms):
            if(molecule.topology.bond_matrix[int(i*(i-1)/2 + atom)]):
                rotate_atom(molecule, i, cos0, sin0, t)
                recursive_fast_rotate1(molecule, i)
                
    init_rotation_axis(molecule, a, b)
    cos0 = np.cos(theta) # cosine of the angle of rotation
    sin0 = np.sin(theta) # sine of the angle of rotation
    t = 1 - cos0
    atom = b-1
    if(a < b):
        recursive_fast_rotate1(molecule, atom)
    #else:
    
def rotate_list(molecule, atom, rlist, k):
    for i in range(atom+1, molecule.num_atoms):
        if(molecule.topology.bond_matrix[int(i*(i-1)/2 + atom)]):
            rlist.append(i)
            rotate_list(molecule, i, rlist, k+1)
            
#def minimize_pot_energy(molecule, atom, d):
#    molecule.x[atom] -= molecule.x[atom]*d
#    molecule.y[atom] -= molecule.y[atom]*d
#    molecule.z[atom] -= molecule.z[atom]*d

def calculate_pot_energy(molecule):
    Vb = 0
    natoms = molecule.num_atoms
    for atom in range(1, natoms):
        for i in range(atom):
            if(molecule.topology.bond_matrix[int(atom*(atom-1)/2 + i)]):
                btype = ('{}-{}').format(molecule.topology.type[atom], 
                        molecule.topology.type[i])  # bond type (e.g. 'c3-c3' or 'c3-hc')
                try:
                    K = molecule.topology.bond_types[btype][0]  # elastic constant of the bond
                    b0 = molecule.topology.bond_types[btype][1] # equilibrium bond length
                except:
                    btype = ('{}-{}').format(molecule.topology.type[i], 
                            molecule.topology.type[atom])
                    try:
                        K = molecule.topology.bond_types[btype][0]
                        b0 = molecule.topology.bond_types[btype][1]
                    except: pass
                b = (((molecule.x[atom] - molecule.x[i])**2) + \
                     ((molecule.y[atom] - molecule.y[i])**2) + \
                     ((molecule.z[atom] - molecule.z[i])**2))**0.5
                d = b-b0
                Vb += K*(d**2) # acumulates the potential energy of the bond
    return Vb

def pot_energy_rotate(molecule, a, b, theta1, theta2):
    init_rotation_axis(molecule, a, b)
    rlist = list()
    degrees = list()
    Vb = []
    if(a < b):
        atom = b-1
        rotate_list(molecule, atom, rlist, 0)
        while(theta1 <= theta2):
            cos0 = np.cos(theta1) # cosine of the angle of rotation
            sin0 = np.sin(theta1) # sine of the angle of rotation
            t = 1 - cos0
            for i in range(len(rlist)):
                rotate_atom(molecule, rlist[i], cos0, sin0, t)
            Vb.append(calculate_pot_energy(molecule))
            degrees.append(theta1)
            theta1 += np.pi/36
        plt.plot(degrees, Vb)
        plt.xlabel('Degrees')
        plt.ylabel('Elastic potential')
        plt.title('Rotation degrees (rad) X Elastic potential (Vb) graphic')
        plt.savefig('elastic_potential_graph_butane_opt.pdf')
        plt.show()
            
molecule = Molecule()

#pdb, psf, out_pdb, out_psf = molecule.get_cmd_line()

molecule.read_psf("butane.psf")
molecule.read_pdb("butane_opt.pdb")
molecule.read_frcmod("butane.frcmod")
#bond_list_C2 = molecule.get_bond_list(2)
#angle_list_C3 = molecule.get_angle_list(3)
#dihed_list_C3 = molecule.get_dihedral_list(3)
#fast_rotate(molecule, 2, 3, np.pi/2)
#Vb = calculate_pot_energy(molecule)
#print("Vb =", Vb)
pot_energy_rotate(molecule, 2, 3, 0, 2*np.pi)
molecule.write_psf("butane1.psf")
molecule.write_pdb("butane_opt_rotated1.pdb")

#for i in range(0, len(bond_list_C3), 2):
#	print(bond_list_C3[i:i+2])
#print("\n")
#for i in range(0, len(angle_list_C3), 3):
#	print(angle_list_C3[i:i+3])
#print("\n")
#for i in range(0, len(dihed_list_C3), 4):
#	print(dihed_list_C3[i:i+4])
#print("\n")
#for i in range(0, len(molecule.topology.bonds), 2):
#	print(molecule.topology.bond_list[i:i+2])
#print("\n")
