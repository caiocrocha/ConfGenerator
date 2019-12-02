import numpy as np
from class_molecule import Molecule

def recursive_rotate_list(molecule, r, c, start, b):
    if(c < len(r) < molecule.num_atoms):
        for i in range(start, molecule.topology.num_bonds*2, 2):
            p = molecule.topology.bond_list[i]
            q = molecule.topology.bond_list[i+1]
            e = r[c]
            if(e == p and q != b):
                r.append(q)
                start = i+2
                break
            elif(e == q and p != b): 
                r.append(p)
                start = i+2
                break
        for i in range(start, molecule.topology.num_bonds*2, 2):
            p = molecule.topology.bond_list[i]
            q = molecule.topology.bond_list[i+1]
            e = r[c]
            if(e == p and q != b):
                r.append(q)
            elif(e == q and p != b): 
                r.append(p)
        recursive_rotate_list(molecule, r, c+1, start, b)

def rotate_list(molecule, r, a, b):
    for i in range(0, molecule.topology.num_bonds*2, 2):
        p = molecule.topology.bond_list[i]
        q = molecule.topology.bond_list[i+1]
        if(p == b and q != a):
            r.append(q)
            start = i+2
            break
        elif(q == b and p != a): 
            r.append(p)
            start = i+2
            break
    for i in range(start, molecule.topology.num_bonds*2, 2):
        p = molecule.topology.bond_list[i]
        q = molecule.topology.bond_list[i+1]
        if(p == b and q != a):
            r.append(q)
        elif(q == b and p != a): 
            r.append(p)
    recursive_rotate_list(molecule, r, 0, start, b)
       
def rotate(molecule, a, b, theta):   # bond a-b, rotate chain from b // theta in radians
    r = list()  # list of atoms to be rotated
    rotate_list(molecule, r, a, b)
    r = list(set(r)) # removes duplicates
    
#    For the construction of the rotation matrix around the axis defined by the bond that
#    connects atoms a and b (vector ab, a to b), there needs to be determined ab's unit vector, 
#    which is mathematicaly defined by its original components divided by the square root
#    of its magnitude.
    
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
    cos0 = np.cos(theta) # cosine of the angle of rotation
    sin0 = np.sin(theta) # sine of the angle of rotation
    t = 1 - cos0
    for i in range(len(r)):
        j = r[i]-1
        x1 = molecule.x[j] - Ax
        y1 = molecule.y[j] - Ay
        z1 = molecule.z[j] - Az
        molecule.x[j] = round(((cos0 + t*(ux**2))*x1 + (t*ux*uy - sin0*uz)*y1 + (t*ux*uz + sin0*uy)*z1), 3) + Ax
        molecule.y[j] = round(((t*ux*uy + sin0*uz)*x1 + (cos0 + t*(uy**2))*y1 + (t*uy*uz - sin0*ux)*z1), 3) + Ay
        molecule.z[j] = round(((t*ux*uz - sin0*uy)*x1 + (t*uy*uz + sin0*ux)*y1 + (cos0 + t*(uz**2))*z1), 3) + Az
        # rotation matrix applied to point
        
def calculate_pot_energy(molecule):
    Vb = 0  # total elastic potential energy in the molecule
    for i in range(0, molecule.topology.num_bonds*2, 2):
        for j in range(len(molecule.topology.bond_types)):
            bond = (molecule.topology.bond_types[j].split('-')) # bond type (e.g. 'c3-c3' or 'c3-hc')
            ind_a = molecule.topology.bond_list[i]-1 # index of atom a
            ind_b = molecule.topology.bond_list[i+1]-1 # index of atom b
            atom_a = molecule.topology.type[ind_a] # atom type (e.g. 'c3' or 'hc')
            atom_b = molecule.topology.type[ind_b]
            if((atom_a == bond[0] and atom_b == bond[1]) or (atom_a == bond[1] and atom_b == bond[0])):
                K = molecule.topology.k_elastic[j] # elastic constant of the bond
                b = (((molecule.x[ind_a]-molecule.x[ind_b])**2)+((molecule.y[ind_a]-molecule.y[ind_b])**2)+((molecule.z[ind_a]-molecule.z[ind_b])**2))**0.5
                # b is the distance between atoms a and b
                b0 = molecule.topology.bond_lengths[j] # b0 is the equilibrium bond length
                Vb += K*((b-b0)**2) # acumulates the potential energy of the bond
                break
    return Vb
        
molecule = Molecule()

#pdb, psf, out_pdb, out_psf = molecule.get_cmd_line()

molecule.read_psf("butane.psf")
molecule.read_pdb("butane_opt.pdb")
#molecule.read_frcmod("butane.frcmod")
bond_list_C3 = molecule.get_bond_list(4)
angle_list_C3 = molecule.get_angle_list(4)
dihed_list_C3 = molecule.get_dihedral_list(4)
#rotate(molecule, 2, 3, np.pi/2)
#Vb = calculate_pot_energy(molecule)
#print("Vb =", Vb)
#molecule.write_psf("butane1.psf")
#molecule.write_pdb("butane1.pdb")

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
