import numpy as np
import math
import matplotlib.pyplot as plt
from class_molecule import Molecule
from rotate import get_rotation_list, init_rotation_axis, rotate_atom

def is_main():
    return __name__ == '__main__'

def angle_Ep_rotate(molecule, a, b, theta, ntimes, pdf_name, 
                      write_pdb=False, pdb_name=None):
    times = ntimes + 1
    degrees = np.zeros(times, dtype='float16')
    Va = np.zeros(times, dtype='float32')
    degrees[0] = 0
    Va[0] = angle_potential(molecule)
    init_rotation_axis(molecule, a, b, theta)
    rlist = get_rotation_list(molecule, a, b)
    if not write_pdb:
        for i in range(1, times):
            for atom in rlist:
                rotate_atom(molecule, atom)
            degrees[i] = theta*i
            Va[i] = angle_potential(molecule)
    else:
        molecule.write_pdb(pdb_name, 'w+', 0)
        for i in range(1, times):
            for atom in rlist:
                rotate_atom(molecule, atom)
            degrees[i] = theta*i
            Va[i] = angle_potential(molecule)
            molecule.write_pdb(pdb_name, 'a', i)
    plt.plot(degrees, Va)
    plt.xlabel('Degrees (rad)')
    plt.ylabel(r'Elastic potential (kcal $mol^{-1} \AA^{-2}$)')
    plt.grid()
    plt.suptitle('Rotation degrees X Elastic potential (Va) graphic')
    plt.title(r'{}$^\circ$ x {} rotations'.format(theta*180/math.pi, ntimes), fontsize=10, loc='right')
    plt.savefig(pdf_name)
    plt.show()
    return Va

def minimize_pot_energy(molecule, fx, fy, fz):
    ''' completa com forca resultante de ligacoes, angulos e diedros '''
    for i in range(molecule.num_atoms):
        fx1 = fx[i]
        fy1 = fy[i]
        fz1 = fz[i]
        N = (fx1*fx1 + fy1*fy1 + fz1*fz1)**0.5 # normal force
        molecule.x[i] += 0.01*fx1/N
        molecule.y[i] += 0.01*fy1/N
        molecule.z[i] += 0.01*fz1/N

def bond_variables(molecule, atom1, atom2):
    bx = molecule.x[atom2] - molecule.x[atom1]
    by = molecule.y[atom2] - molecule.y[atom1]
    bz = molecule.z[atom2] - molecule.z[atom1]
    b = (bx*bx + by*by + bz*bz)**0.5
    ux = bx/b
    uy = by/b
    uz = bz/b
    btype = ('{}-{}'.format(molecule.topology.type[atom1], 
                    molecule.topology.type[atom2])
             )  # bond type (e.g. 'c3-c3' or 'c3-hc')
    try:
        Kb = molecule.topology.bond_types[btype][0] # elastic constant of the bond
        b0 = molecule.topology.bond_types[btype][1] # equilibrium bond length
    except:
        btype = ('{}-{}'.format(molecule.topology.type[atom2], 
                molecule.topology.type[atom1])
                 )
        try:
            Kb = molecule.topology.bond_types[btype][0]
            b0 = molecule.topology.bond_types[btype][1]
        except: pass
    d = b-b0
    # d: difference between the present distance between the two atoms 
    # and the equilibrium distance
    return ux, uy, uz, b, b0, d, Kb

def bond_potential(molecule):
    nbonds = molecule.topology.num_bonds
    Vb = np.zeros(nbonds, dtype='float32')
    count = 0
    for i in range(0, nbonds*2, 2):
        ux, uy, uz, b, b0, d, Kb = bond_variables(molecule, 
            atom1 = molecule.topology.bond_list[i]-1, 
            atom2 = molecule.topology.bond_list[i+1]-1)
        Vb[count] = Kb*d*d        
        # print('count =', count)
        # print('{:>10s}{:>10s}{:>10s}'.format('fx', 'fy', 'fz'))
        # print('{:2d}:{:7.3f}{:10.3f}{:10.3f}'.format(atom1+1, fx1, fy1, fz1))
        # print('{:2d}:{:7.3f}{:10.3f}{:10.3f}'.format(atom2+1, fx2, fy2, fz2))
        # for i in range(30):
        #     print('-', end='')
        # print('\n{:2d}:{:7.3f}{:10.3f}{:10.3f}'.format(atom1+1, fx[atom1], fy[atom1], fz[atom1]))
        # print('{:2d}:{:7.3f}{:10.3f}{:10.3f}'.format(atom2+1, fx[atom2], fy[atom2], fz[atom2]))
        # print('END', end='')
        # for i in range(27):
        #     print('_', end='')
        # print('\n')
        count += 1
    return Vb

def bond_force(molecule, fx, fy, fz):
    for i in range(0, molecule.topology.num_bonds*2, 2):
        atom1 = molecule.topology.bond_list[i]-1
        atom2 = molecule.topology.bond_list[i+1]-1
        ux, uy, uz, b, b0, d, Kb = bond_variables(molecule, atom1, atom2)
        force = -2*Kb*d  # first derivative of the bond potential Vb
        fx1 = force*ux  # decomposition of the force on the X axis
        fy1 = force*uy  # decomposition of the force on the Y axis
        fz1 = force*uz  # decomposition of the force on the Z axis
        ''' ta como bx, by, bz; nao seria ux, uy, uz ?? '''
        fx2 = -fx1
        fy2 = -fy1
        fz2 = -fz1
        fx[atom1] += fx1
        fy[atom1] += fy1
        fz[atom1] += fz1
        fx[atom2] += fx2
        fy[atom2] += fy2
        fz[atom2] += fz2
        
def cross_product_3d(molecule, v1x, v2x, v1y, v2y, v1z, v2z):
    # The cross product is the determinant of the matrix below:
    #       |  i  j  k |
    #   det | x1 y1 z1 |
    #       | x2 y2 z2 |
    #   = ( y1 * z2 - z1 * y2 ) * i
    #   + ( z1 * x2 - x1 * z2 ) * j
    #   + ( x1 * y2 - y1 * x2 ) * k
    v3x = v1y*v2z - v1z*v2y
    v3y = v1z*v2x - v1x*v2z
    v3z = v1x*v2y - v1y*v2x
    return v3x, v3y, v3z

def angle_pot_variables(molecule, atom1, atom2, atom3):
    # global x1, y1, z1, x2, y2, z2, x3, y3, z3
    # global v21x, v21y, v21z, v23x, v23y, v23z, n1, n2
    
    x1 = molecule.x[atom1]
    y1 = molecule.y[atom1]
    z1 = molecule.z[atom1]
    x2 = molecule.x[atom2]
    y2 = molecule.y[atom2]
    z2 = molecule.z[atom2]
    x3 = molecule.x[atom3]
    y3 = molecule.y[atom3]
    z3 = molecule.z[atom3]
    
    v21x = x1 - x2
    v23x = x3 - x2
    v21y = y1 - y2
    v23y = y3 - y2
    v21z = z1 - z2
    v23z = z3 - z2
    norm1 = (v21x*v21x + v21y*v21y + v21z*v21z)**0.5 # norm of vector v21 (or v12)
    norm2 = (v23x*v23x + v23y*v23y + v23z*v23z)**0.5 # norm of vector v23 (or v32)
    
    atype = ('{}-{}-{}'.format(molecule.topology.type[atom1], 
                    molecule.topology.type[atom2], molecule.topology.type[atom3])
             )  # angle type (e.g. 'c3-c3-hc' or 'c3-c3-c3')
    try:
        Ka = molecule.topology.angle_types[atype][0]    # angle elastic constant
        a0 = molecule.topology.angle_types[atype][1]    # equilibrium angle
    except:
        atype = ('{}-{}-{}'.format(molecule.topology.type[atom3], 
                    molecule.topology.type[atom2], molecule.topology.type[atom1])
                 )
        try:
            Ka = molecule.topology.angle_types[atype][0]
            a0 = molecule.topology.angle_types[atype][1]
        except: pass
    
    angle = math.acos((v21x*v23x + v21y*v23y + v21z*v23z)/(norm1*norm2))
    delta_angle = angle - a0
    return delta_angle, Ka

def angle_force_variables(molecule, atom1, atom2, atom3):
    x1 = molecule.x[atom1]
    y1 = molecule.y[atom1]
    z1 = molecule.z[atom1]
    x2 = molecule.x[atom2]
    y2 = molecule.y[atom2]
    z2 = molecule.z[atom2]
    x3 = molecule.x[atom3]
    y3 = molecule.y[atom3]
    z3 = molecule.z[atom3]
    
    v21x = x1 - x2
    v23x = x3 - x2
    v21y = y1 - y2
    v23y = y3 - y2
    v21z = z1 - z2
    v23z = z3 - z2
    norm1 = (v21x*v21x + v21y*v21y + v21z*v21z)**0.5 # norm of vector v21 (or v12)
    norm2 = (v23x*v23x + v23y*v23y + v23z*v23z)**0.5 # norm of vector v23 (or v32)
    
    atype = ('{}-{}-{}'.format(molecule.topology.type[atom1], 
                    molecule.topology.type[atom2], molecule.topology.type[atom3])
             )  # angle type (e.g. 'c3-c3-hc' or 'c3-c3-c3')
    try:
        Ka = molecule.topology.angle_types[atype][0]    # angle elastic constant
        a0 = molecule.topology.angle_types[atype][1]    # equilibrium angle
    except:
        atype = ('{}-{}-{}'.format(molecule.topology.type[atom3], 
                    molecule.topology.type[atom2], molecule.topology.type[atom1])
                 )
        try:
            Ka = molecule.topology.angle_types[atype][0]
            a0 = molecule.topology.angle_types[atype][1]
        except: pass
    
    angle = math.acos((v21x*v23x + v21y*v23y + v21z*v23z)/(norm1*norm2))
    delta_angle = angle - a0
    
    nx, ny, nz = cross_product_3d(
        molecule, v21x, v23x, v21y, v23y, v21z, v23z)
    # normal vector of the plane determined by vectors v21 and v23
    vR1x, vR1y, vR1z = cross_product_3d(
        molecule, v21x, nx, v21y, ny, v21z, nz)
    # vR1: resulting vector (orthogonal to v21)
    vR2x, vR2y, vR2z = cross_product_3d(
        molecule, v23x, nx, v23y, ny, v23z, nz)
    # vR2: resulting vector (orthogonal to v23)
    M1 = (vR1x*vR1x + vR1y*vR1y + vR1z*vR1z)**0.5
    M2 = (vR2x*vR2x + vR2y*vR2y + vR2z*vR2z)**0.5
    # the following are the normalization of vR1 and vR2:
    p1x = vR1x/M1
    p1y = vR1y/M1
    p1z = vR1z/M1
    p2x = vR2x/M2
    p2y = vR2y/M2
    p2z = vR2z/M2
    return p1x, p1y, p1z, p2x, p2y, p2z, delta_angle, Ka
        
def angle_potential(molecule):
    nangles = molecule.topology.num_angles
    Va = np.zeros(nangles, dtype='float32')
    count = 0
    for i in range(0, nangles*3, 3):
        delta_angle, Ka = angle_pot_variables(molecule, 
            atom1 = molecule.topology.angle_list[i]-1, 
            atom2 = molecule.topology.angle_list[i+1]-1, 
            atom3 = molecule.topology.angle_list[i+2]-1)
        Va[count] = Ka*delta_angle*delta_angle
        count += 1
    return Va
        
def angle_force(molecule, fx, fy, fz):
    for i in range(0, molecule.topology.num_angles*3, 3):
        atom1 = molecule.topology.angle_list[i]-1
        atom2 = molecule.topology.angle_list[i+1]-1
        atom3 = molecule.topology.angle_list[i+2]-1
        p1x, p1y, p1z, p2x, p2y, p2z, norm1, norm2, delta_angle, Ka = angle_force_variables(
            molecule, atom1, atom2, atom3)
        force = -2*Ka*delta_angle    # first derivative of the angle potential Va
        u1 = force/norm1
        # force multiplied by the partial derivative of Va according to position (x1, y1, z1)
        u3 = force/norm2
        # force multiplied by the partial derivative of Va according to position (x3, y3, z3)
        fx1 = p1x*u1
        fy1 = p1y*u1
        fz1 = p1z*u1
        fx2 = p2x*u3
        fy2 = p2y*u3
        fz2 = p2z*u3
        fx[atom1] += fx1
        fy[atom1] += fy1
        fz[atom1] += fz1
        fx[atom3] += fx2
        fy[atom3] += fy2
        fz[atom3] += fz2
        fx[atom2] += -fx1 - fx2
        ''' fx1 ou fx[atom1] ?? '''
        fy[atom2] += -fy1 - fy2
        fz[atom2] += -fz1 - fz2
        
def dihedral_pot_variables(molecule, atom1, atom2, atom3, atom4):
    x1 = molecule.x[atom1]
    y1 = molecule.y[atom1]
    z1 = molecule.z[atom1]
    x2 = molecule.x[atom2]
    y2 = molecule.y[atom2]
    z2 = molecule.z[atom2]
    x3 = molecule.x[atom3]
    y3 = molecule.y[atom3]
    z3 = molecule.z[atom3]
    x4 = molecule.x[atom4]
    y4 = molecule.y[atom4]
    z4 = molecule.z[atom4]
    
    v21x = x1 - x2
    v21y = y1 - y2
    v21z = z1 - z2
    v23x = x3 - x2
    v23y = y3 - y2
    v23z = z3 - z2
    v34x = x4 - x3
    v34y = y4 - y3
    v34z = z4 - z3
    norm1 = (v21x*v21x + v21y*v21y + v21z*v21z)**0.5 # norm of vector v21 (or v12)
    norm2 = (v23x*v23x + v23y*v23y + v23z*v23z)**0.5 # norm of vector v23 (or v32)
    norm3 = (v34x*v34x + v34y*v34y + v34z*v34z)**0.5 # norm of vector v34 (or v43)
    
    n1x, n1y, n1z = cross_product_3d(
        molecule, v21x, v23x, v21y, v23y, v21z, v23z)
    # normal vector of the plane determined by vectors v21 and v23
    n2x, n2y, n2z = cross_product_3d(
        molecule, -v23x, v34x, -v23y, v34y, -v23z, v34z)
    # normal vector of the plane determined by vectors v32 and v34 (v32 = -v23)
    
    M1 = (n1x*n1x + n1y*n1y + n1z*n1z)**0.5 # magnitude of n1
    M2 = (n2x*n2x + n2y*n2y + n2z*n2z)**0.5 # magnitude of n2
    dihed_angle = math.acos((n1x*n2x + n1y*n2y + n1z*n2z)/(M1*M2))
    
    return dihed_angle
        
def dihedral_force_variables(molecule, atom1, atom2, atom3, atom4):
    x1 = molecule.x[atom1]
    y1 = molecule.y[atom1]
    z1 = molecule.z[atom1]
    x2 = molecule.x[atom2]
    y2 = molecule.y[atom2]
    z2 = molecule.z[atom2]
    x3 = molecule.x[atom3]
    y3 = molecule.y[atom3]
    z3 = molecule.z[atom3]
    x4 = molecule.x[atom4]
    y4 = molecule.y[atom4]
    z4 = molecule.z[atom4]
    
    v21x = x1 - x2
    v21y = y1 - y2
    v21z = z1 - z2
    v23x = x3 - x2
    v23y = y3 - y2
    v23z = z3 - z2
    v34x = x4 - x3
    v34y = y4 - y3
    v34z = z4 - z3
    norm1 = (v21x*v21x + v21y*v21y + v21z*v21z)**0.5 # norm of vector v21 (or v12)
    norm2 = (v23x*v23x + v23y*v23y + v23z*v23z)**0.5 # norm of vector v23 (or v32)
    norm3 = (v34x*v34x + v34y*v34y + v34z*v34z)**0.5 # norm of vector v34 (or v43)
    
    n1x, n1y, n1z = cross_product_3d(
        molecule, v21x, v23x, v21y, v23y, v21z, v23z)
    # normal vector of the plane determined by vectors v21 and v23
    n2x, n2y, n2z = cross_product_3d(
        molecule, -v23x, v34x, -v23y, v34y, -v23z, v34z)
    # normal vector of the plane determined by vectors v32 and v34 (v32 = -v23)
    
    M1 = (n1x*n1x + n1y*n1y + n1z*n1z)**0.5
    M2 = (n2x*n2x + n2y*n2y + n2z*n2z)**0.5
    dihed_angle = math.acos((n1x*n2x + n1y*n2y + n1z*n2z)/(M1*M2))
    
    p1x = n1x/M1
    p1y = n1y/M1
    p1z = n1z/M1
    p2x = n2x/M2
    p2y = n2y/M2
    p2z = n2z/M2
    # p1 and p2 are the normalization of n1 and n2
    
    force = 0.5*(A1*math.sin(dihed_angle) - 2*A2*math.sin(2*dihed_angle) + 
                  3*A3*math.sin(3*dihed_angle)) # first derivative of the dihedral potential Vd
    u1 = force/(norm1*math.sin(theta1))
    # force multiplied by the partial derivative of theta1 according to position (x1, y1, z1)
    u4 = force/(norm3*math.sin(theta2))
    # force multiplied by the partial derivative of theta2 according to position (x4, y4, z4)
    f1x = p1x*u1
    f1y = p1y*u1
    f1z = p1z*u1
    f4x = p2x*u4
    f4y = p2y*u4
    f4z = p2z*u4
    
    ox = (x2 + x3)/2
    oy = (y2 + y3)/2
    oz = (z2 + z3)/2
    # The point o is the center of bond v23
    vo3x = x3 - ox
    vo3y = y3 - oy
    vo3z = z3 - oz
    # vector parting from o to (x3, y3, z3)
    
    const = -1/(vo3x*vo3x + vo3y*vo3y + vo3z*vo3z)   # inverse square of the norm of vector vo3
    n_vo3_f4x, n_vo3_f4y, n_vo3_f4z = cross_product_3d(molecule, vo3x, f4x, vo3y, f4y, vo3z, f4z)
    n_v34_f4x, n_v34_f4y, n_v34_f4z = cross_product_3d(molecule, v34x, f4x, v34y, f4y, v34z, f4z)
    n_v21_f4x, n_v21_f4y, n_v21_f4z = cross_product_3d(molecule, v21x, f4x, v21y, f4y, v21z, f4z)
    t3x = const*(n_vo3_f4x + 0.5*n_v34_f4x + 0.5*n_v21_f4x)
    t3y = const*(n_vo3_f4y - 0.5*n_v34_f4y - 0.5*n_v21_f4y)
    t3z = const*(n_vo3_f4z - 0.5*n_v34_f4z - 0.5*n_v21_f4z)
    # t3 is the cross product (vo3 by f3)
    
    f3x, f3y, f3z = cross_product_3d(molecule, t3x, vo3x, t3y, vo3y, t3z, vo3z)
    f2x = -f1x - f3x - f4x
    f2y = -f1y - f3y - f4y
    f2z = -f1z - f3z - f4z
    return f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z, f4x, f4y, f4z

def dihedral_potential(molecule):
    ndihed = molecule.topology.num_dihedrals
    Vd = np.zeros(ndihed, dtype='float32')
    count = 0
    # Vd = 0
    for i in range(0, ndihed*4, 4):
        dihed_angle = dihedral_pot_variables(molecule, 
            atom1 = molecule.topology.dihedral_list[i]-1, 
            atom2 = molecule.topology.dihedral_list[i+1]-1, 
            atom3 = molecule.topology.dihedral_list[i+2]-1, 
            atom4 = molecule.topology.dihedral_list[i+3]-1)
        Vd[count] = 0.5*(
            A1*(1 + math.cos(dihed_angle)) + A2*(1 - math.cos(2*dihed_angle)) + 
            A3*(1 + math.cos(3*dihed_angle)) + A4)
        count += 1
        # Vd += Ka*delta_angle*delta_angle
    return Vd

def dihedral_force(molecule, fx, fy, fz):
    for i in range(0, molecule.topology.num_dihedrals*4, 4):
        atom1 = molecule.topology.dihedral_list[i]-1
        atom2 = molecule.topology.dihedral_list[i+1]-1
        atom3 = molecule.topology.dihedral_list[i+2]-1
        atom4 = molecule.topology.dihedral_list[i+3]-1
        fx1, fy1, fz1, fx2, fy2, fz2, fx3, fy3, fz3, fx4, fy4, fz4 = dihedral_force_variables(
            molecule, atom1, atom2, atom3, atom4)
        fx[atom1] += fx1
        fy[atom1] += fy1
        fz[atom1] += fz1
        fx[atom2] += fx2
        fy[atom2] += fy2
        fz[atom2] += fz2
        fx[atom3] += fx3
        fy[atom3] += fy3
        fz[atom3] += fz3
        fx[atom4] += fx4
        fy[atom4] += fy4
        fz[atom4] += fz4

##############################################################################
if is_main():
    
    molecule = Molecule()
    
    # pdb, psf, out_pdb, out_psf = molecule.get_cmd_line()
    
    molecule.read_psf('butane.psf')
    molecule.read_pdb('sqm.pdb')
    molecule.read_frcmod('butane.frcmod')
    # Va = angle_Ep_rotate(molecule, a=2, b=3, theta=np.pi/18, ntimes=36, 
    #                        pdf_name='Ep_sqm.pdf', write_pdb=False, pdb_name='sqm_rotated.pdb')
    Va = angle_potential(molecule)
    # fx = np.zeros(molecule.num_atoms, dtype='float32')
    # fy = np.zeros(molecule.num_atoms, dtype='float32')
    # fz = np.zeros(molecule.num_atoms, dtype='float32')
    # bond_force(molecule, fx, fy, fz)
    # angle_force(molecule, fx, fy, fz)