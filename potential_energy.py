import numpy as np
import math
import matplotlib.pyplot as plt
from class_molecule import Molecule
from rotate import get_rotation_list, init_rotation_axis, set_rotation_angle, rotate_atom
from cross_product_3d import cross_product_3d
from plot_Ep import plot_Ep

def is_main():
    return __name__ == '__main__'

# %%

def heuristic_rotate(molecule, rlist, theta, ntimes, arquivo, atom1, atom2, atom3, atom4, vez):
    init_rotation_axis(molecule, atom2 + 1, atom3 + 1, theta)
    lowest, dihed0 = _angle_dihedral_potential_gaff(molecule, atom1, atom2, atom3, atom4)
    dihed0 *= 180 / np.pi
    aux1 = dihedral_potential(molecule)
    comeback = 0
    for i in range(ntimes - 1):
        for atom in rlist:
            rotate_atom(molecule, atom)
        # molecule.write_pdb(arquivo, 'a', vez)
        vez += 1
        Vd, dihed1 = _angle_dihedral_potential_gaff(molecule, atom1, atom2, atom3, atom4)
        dihed1 *= 180/np.pi
        aux2 = dihedral_potential(molecule)
        if aux2 < aux1:
            lowest = Vd
            comeback = i + 1
            aux1 = aux2
    set_rotation_angle(theta * (1 + comeback))
    for atom in rlist:
        rotate_atom(molecule, atom)
    # molecule.write_pdb(arquivo, 'a', vez)
    vez += 1

def heuristic_conformation(molecule, theta, arquivo):
    ntimes = int(2 * np.pi / theta)
    global vez
    vez = 1
    # molecule.write_pdb(arquivo, 'w+', 0)
    for i in range(0, molecule.topology.num_dihedrals, 4):
        atom1 = molecule.topology.dihedral_list[i] - 1
        atom2 = molecule.topology.dihedral_list[i + 1] - 1
        atom3 = molecule.topology.dihedral_list[i + 2] - 1
        atom4 = molecule.topology.dihedral_list[i + 3] - 1
        if molecule.topology.type[atom1] == 'hc' or molecule.topology.type[atom4] == 'hc':
            continue
        rlist = get_rotation_list(molecule, atom2 + 1, atom3 + 1)
        heuristic_rotate(molecule, rlist, theta, ntimes, arquivo, atom1, atom2, atom3, atom4, vez)
    molecule.write_mol2(arquivo)

def dihedral_Ep_rotate(molecule, a, b, theta, ntimes, pdf_name,
                       write_pdb=False, pdb_name=None):
    times = ntimes + 1
    Vd = np.zeros(times, dtype='float32')
    dihed_angle = np.zeros(times, dtype='float32')
    Vd[0], dihed_angle[0] = _dihedral_potential(molecule)
    plt.scatter(dihed_angle[0], Vd[0], marker='.', color='royalblue')
    init_rotation_axis(molecule, a, b, theta)
    rlist = get_rotation_list(molecule, a, b)
    if not write_pdb:
        for i in range(1, times):
            for atom in rlist:
                rotate_atom(molecule, atom)
            Vd[i], dihed_angle[i] = _dihedral_potential(molecule)
            dihed_angle[i] = dihed_angle[i]*180/np.pi + 180
            plt.scatter(dihed_angle[i], Vd[i], marker='.', color='royalblue')
    else:
        molecule.write_pdb(pdb_name, 'w+', 0)
        for i in range(1, times):
            for atom in rlist:
                rotate_atom(molecule, atom)
            Vd[i], dihed_angle[i] = _dihedral_potential(molecule)
            dihed_angle[i] = dihed_angle[i]*180/np.pi + 180
            plt.scatter(dihed_angle[i], Vd[i], marker='.', color='royalblue')
            molecule.write_pdb(pdb_name, 'a', i)
    plot_Ep(theta, ntimes, pdf_name)
    return dihed_angle, Vd

def minimize_pot_energy(molecule, fx, fy, fz):
    for i in range(molecule.num_atoms):
        fx1 = fx[i]
        fy1 = fy[i]
        fz1 = fz[i]
        N = (fx1 * fx1 + fy1 * fy1 + fz1 * fz1) ** 0.5  # normal force
        molecule.x[i] -= 0.01 * fx1 / N
        molecule.y[i] -= 0.01 * fy1 / N
        molecule.z[i] -= 0.01 * fz1 / N

def bond_pot_variables(molecule, atom1, atom2):
    bx = molecule.x[atom2] - molecule.x[atom1]
    by = molecule.y[atom2] - molecule.y[atom1]
    bz = molecule.z[atom2] - molecule.z[atom1]
    b = (bx * bx + by * by + bz * bz) ** 0.5  # distance between atoms 1 and 2
    btype = ('{}-{}'.format(molecule.topology.type[atom1],
                            molecule.topology.type[atom2])
             )  # bond type (e.g. 'c3-c3' or 'c3-hc')
    try:
        Kb = molecule.topology.bond_types[btype][0]  # elastic constant of the bond
        b0 = molecule.topology.bond_types[btype][1]  # equilibrium bond length
    except:
        btype = ('{}-{}'.format(molecule.topology.type[atom2],
                                molecule.topology.type[atom1])
                 )
        try:
            Kb = molecule.topology.bond_types[btype][0]
            b0 = molecule.topology.bond_types[btype][1]
        except:
            pass
    d = b - b0
    # d: difference between the present distance between the two atoms
    # and the equilibrium distance
    return bx, by, bz, b, b0, Kb, d

def bond_force_variables(molecule, atom1, atom2):
    bx, by, bz, b, b0, Kb, d = bond_pot_variables(molecule, atom1, atom2)
    if b != 0:
        ux = bx / b
        uy = by / b
        uz = bz / b
    else:
        ux = 0
        uy = 0
        uz = 0
    return bx, by, bz, b, b0, Kb, d, ux, uy, uz

def bond_potential(molecule):
    nbonds = molecule.topology.num_bonds
    Vb = 0
    # Vb = np.zeros(nbonds, dtype='float32')
    # count = 0
    for i in range(0, nbonds * 2, 2):
        bx, by, bz, b, b0, Kb, d = bond_pot_variables(molecule,
                                                      atom1=molecule.topology.bond_list[i] - 1,
                                                      atom2=molecule.topology.bond_list[i + 1] - 1)
        Vb += Kb * d * d
        # Vb[count] = Kb*d*d
        # count += 1
    return Vb

def bond_force(molecule, fx, fy, fz):
    for i in range(0, molecule.topology.num_bonds * 2, 2):
        atom1 = molecule.topology.bond_list[i] - 1
        atom2 = molecule.topology.bond_list[i + 1] - 1
        bx, by, bz, b, b0, Kb, d, ux, uy, uz = bond_force_variables(molecule, atom1, atom2)
        force = -2 * Kb * d  # first derivative of the bond potential Vb
        fx1 = force * ux  # decomposition of the force on the X axis
        fy1 = force * uy  # decomposition of the force on the Y axis
        fz1 = force * uz  # decomposition of the force on the Z axis
        fx2 = -fx1
        fy2 = -fy1
        fz2 = -fz1
        fx[atom1] += fx1
        fy[atom1] += fy1
        fz[atom1] += fz1
        fx[atom2] += fx2
        fy[atom2] += fy2
        fz[atom2] += fz2

def _angle_pot_variables(molecule, atom1, atom2, atom3):
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
    norm1 = (v21x * v21x + v21y * v21y + v21z * v21z) ** 0.5  # norm of vector v21 (or v12)
    norm2 = (v23x * v23x + v23y * v23y + v23z * v23z) ** 0.5  # norm of vector v23 (or v32)

    atype = ('{}-{}-{}'.format(molecule.topology.type[atom1],
                               molecule.topology.type[atom2], molecule.topology.type[atom3])
             )  # angle type (e.g. 'c3-c3-hc' or 'c3-c3-c3')
    try:
        Ka = molecule.topology.angle_types[atype][0]  # angle elastic constant
        a0 = molecule.topology.angle_types[atype][1]  # equilibrium angle
    except:
        atype = ('{}-{}-{}'.format(molecule.topology.type[atom3],
                                   molecule.topology.type[atom2], molecule.topology.type[atom1])
                 )
        try:
            Ka = molecule.topology.angle_types[atype][0]
            a0 = molecule.topology.angle_types[atype][1]
        except:
            raise AttributeError('Could not find corresponding angle type')
    return a0, Ka, v21x, v21y, v21z, v23x, v23y, v23z, norm1, norm2

def angle_pot_variables(molecule, atom1, atom2, atom3):
    a0, Ka, v21x, v21y, v21z, v23x, v23y, v23z, norm1, norm2 = _angle_pot_variables(molecule, atom1, atom2, atom3)
    if norm1 != 0 and norm2 != 0:
        angle = math.acos((v21x * v23x + v21y * v23y + v21z * v23z) / (norm1 * norm2))
    else:
        angle = 0
    return angle - a0, Ka

def angle_force_variables(molecule, atom1, atom2, atom3):
    a0, Ka, v21x, v23x, v21y, v23y, v21z, v23z, norm1, norm2 = _angle_pot_variables(
        molecule, atom1, atom2, atom3)
    nx, ny, nz = cross_product_3d(
        molecule, v21x, v23x, v21y, v23y, v21z, v23z)
    # normal vector of the plane determined by vectors v21 and v23
    vR1x, vR1y, vR1z = cross_product_3d(
        molecule, v21x, nx, v21y, ny, v21z, nz)
    # vR1: resulting vector (orthogonal to v21)
    vR2x, vR2y, vR2z = cross_product_3d(
        molecule, v23x, nx, v23y, ny, v23z, nz)
    # vR2: resulting vector (orthogonal to v23)
    M1 = (vR1x * vR1x + vR1y * vR1y + vR1z * vR1z) ** 0.5
    M2 = (vR2x * vR2x + vR2y * vR2y + vR2z * vR2z) ** 0.5
    # the following are the normalization of vR1 and vR2:
    if M1 != 0:
        p1x = vR1x / M1
        p1y = vR1y / M1
        p1z = vR1z / M1
    else:
        p1x = 0
        p1y = 0
        p1z = 0
    if M2 != 0:
        p2x = vR2x / M2
        p2y = vR2y / M2
        p2z = vR2z / M2
    else:
        p2x = 0
        p2y = 0
        p2z = 0
    return a0, Ka, v21x, v23x, v21y, v23y, v21z, v23z, norm1, norm2, p1x, p1y, p1z, p2x, p2y, p2z

def angle_potential(molecule):
    nangles = molecule.topology.num_angles
    Va = 0
    for i in range(0, nangles * 3, 3):
        delta_angle, Ka = angle_pot_variables(
            molecule,
            atom1=molecule.topology.angle_list[i] - 1,
            atom2=molecule.topology.angle_list[i + 1] - 1,
            atom3=molecule.topology.angle_list[i + 2] - 1)
        Va += Ka * delta_angle * delta_angle
    return Va

def angle_force(molecule, fx, fy, fz):
    for i in range(0, molecule.topology.num_angles * 3, 3):
        atom1 = molecule.topology.angle_list[i] - 1
        atom2 = molecule.topology.angle_list[i + 1] - 1
        atom3 = molecule.topology.angle_list[i + 2] - 1
        a0, Ka, v21x, v23x, v21y, v23y, v21z, v23z, norm1, norm2, p1x, p1y, p1z, p2x, p2y, p2z = angle_force_variables(
            molecule, atom1, atom2, atom3)
        if norm1 != 0 and norm2 != 0:
            dot = (v21x * v23x + v21y * v23y + v21z * v23z) / (norm1 * norm2)
            if dot > 1:
                dot = 1
            elif dot < -1:
                dot = -1
            angle = math.acos(dot)
            force = -2 * Ka * (angle - a0)  # first derivative of the angle potential Va
        else:
            angle = 0
            force = -2 * Ka * (angle - a0)
        if norm1 != 0:
            u1 = force / norm1
            # force multiplied by the partial derivative of Va according to position (x1, y1, z1)
        else:
            u1 = 0
        if norm2 != 0:
            u3 = force / norm2
            # force multiplied by the partial derivative of Va according to position (x3, y3, z3)
        else:
            u3 = 0
        fx1 = p1x * u1
        fy1 = p1y * u1
        fz1 = p1z * u1
        fx2 = p2x * u3
        fy2 = p2y * u3
        fz2 = p2z * u3
        fx[atom1] += fx1
        fy[atom1] += fy1
        fz[atom1] += fz1
        fx[atom3] += fx2
        fy[atom3] += fy2
        fz[atom3] += fz2
        fx[atom2] += -fx1 - fx2
        fy[atom2] += -fy1 - fy2
        fz[atom2] += -fz1 - fz2

def get_dihed_type(molecule, atom1, atom2, atom3, atom4):
    dtype = ('{}-{}-{}-{}'.format(molecule.topology.type[atom1],
                                  molecule.topology.type[atom2],
                                  molecule.topology.type[atom3],
                                  molecule.topology.type[atom4])
             )  # dihedral type (e.g. 'c3-c3-c3-c3')
    return dtype

def new_dihedral(molecule, atom1, atom2, atom3, atom4):
    from time import time
    # Reference: https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    '''
    time0 = time()
    p0x = molecule.x[atom1]
    p0y = molecule.y[atom1]
    p0z = molecule.z[atom1]
    p1x = molecule.x[atom2]
    p1y = molecule.y[atom2]
    p1z = molecule.z[atom2]
    p2x = molecule.x[atom3]
    p2y = molecule.y[atom3]
    p2z = molecule.z[atom3]
    p3x = molecule.x[atom4]
    p3y = molecule.y[atom4]
    p3z = molecule.z[atom4]

    b0x = p0x - p1x
    b0y = p0y - p1y
    b0z = p0z - p1z
    b1x = p2x - p1x
    b1y = p2y - p1y
    b1z = p2z - p1z
    b2x = p3x - p2x
    b2y = p3y - p2y
    b2z = p3z - p2z

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    normb1 = (b1x*b1x + b1y*b1y + b1z*b1z)**0.5
    b1x /= normb1
    b1y /= normb1
    b1z /= normb1

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    normb0 = (b0x*b0x + b0y*b0y + b0z*b0z)**0.5
    normb2 = (b2x*b2x + b2y*b2y + b2z*b2z)**0.5

    dotv = (b0x * b1x + b0y * b1y + b0z * b1z) / (normb0 * normb1)
    vx = b0x - dotv * b1x
    vy = b0y - dotv * b1y
    vz = b0z - dotv * b1z
    dotw = (b2x * b1x + b2y * b1y + b2z * b1z) / (normb2 * normb1)
    wx = b2x - dotw * b1x
    wy = b2y - dotw * b1y
    wz = b2z - dotw * b1z

    normv = (vx * vx + vy * vy + vz * vz) ** 0.5
    normw = (wx * wx + wy * wy + wz * wz) ** 0.5
    x = (vx*wx + vy*wy + vz*wz)/(normv*normw)

    cx, cy, cz = cross_product_3d(b1x, vx, b1y, vy, b1z, vz)
    normc = (cx * cx + cy * cy + cz * cz) ** 0.5
    y = (cx*wx + cy*wy + cz*wz)/(normc*normw)

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    time1 = time() - time0
    return np.arctan2(y, x)

    '''
    time0 = time()
    p0 = np.array([molecule.x[atom1], molecule.y[atom1], molecule.z[atom1]])
    p1 = np.array([molecule.x[atom2], molecule.y[atom2], molecule.z[atom2]])
    p2 = np.array([molecule.x[atom3], molecule.y[atom3], molecule.z[atom3]])
    p3 = np.array([molecule.x[atom4], molecule.y[atom4], molecule.z[atom4]])

    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    time1 = time() - time0
    return np.arctan2(y, x)

def _dihedral_pot_variables(molecule, atom1, atom2, atom3, atom4, dtype):
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
    v32x = -v23x
    v32y = -v23y
    v32z = -v23z

    try:
        length = len(molecule.topology.dihedral_types[dtype])
        A1 = np.zeros(length)
        A2 = np.zeros(length)
        A3 = np.zeros(length)
        A4 = np.zeros(length)
        for i in range(length):
            A1[i] = molecule.topology.dihedral_types[dtype][i][0]
            A2[i] = molecule.topology.dihedral_types[dtype][i][1]
            A3[i] = molecule.topology.dihedral_types[dtype][i][2]
            A4[i] = molecule.topology.dihedral_types[dtype][i][3]
    except:
        dtype = ('{}-{}-{}-{}'.format(molecule.topology.type[atom4],
                                      molecule.topology.type[atom3],
                                      molecule.topology.type[atom2],
                                      molecule.topology.type[atom1])
                 )  # dihedral type (e.g. 'c3-c3-c3-c3')
        try:
            length = len(molecule.topology.dihedral_types[dtype])
            A1 = np.zeros(length)
            A2 = np.zeros(length)
            A3 = np.zeros(length)
            A4 = np.zeros(length)

            for i in range(length):
                A1[i] = molecule.topology.dihedral_types[dtype][i][0]
                A2[i] = molecule.topology.dihedral_types[dtype][i][1]
                A3[i] = molecule.topology.dihedral_types[dtype][i][2]
                A4[i] = molecule.topology.dihedral_types[dtype][i][3]
        except:
            raise AttributeError('Could not find corresponding dihedral type')

    n1x, n1y, n1z = cross_product_3d(v21x, v23x, v21y, v23y, v21z, v23z)
    # normal vector of the plane determined by vectors v21 and v23
    n2x, n2y, n2z = cross_product_3d(v32x, v34x, v32y, v34y, v32z, v34z)
    # normal vector of the plane determined by vectors v32 and v34
    '''
    nx, ny, nz = cross_product_3d(molecule, n1x, n2x, n1y, n2y, n1z, n2z)
    # normal vector of the plane determined by n1 and n2
    M = (nx*nx + ny*ny + nz*nz)**0.5    # magnitude of normal vector n
    if M != 0:
        det = (nx*abs(nx) + ny*abs(ny) + nz*abs(nz))/M
        # Let u = abs(n)/M (normalization of vector n). The determinant det(n1, n2, u)
        # is proportional to the sine of the angle between vectors n1 and n2. It can be
        # expressed as the triple product between n1, n2 and u. Triple product:
        # (n1 x n2) * u
        dot = n1x * n2x + n1y * n2y + n1z * n2z
        # the dot product between n1 and n2 is proportional to the cosine of the angle
        # dihed_angle = math.atan2(det, dot)
        dihed_angle = math.atan2(det, dot)
    else:
        dihed_angle = 0
    '''
    dihed_angle = new_dihedral(molecule, atom1, atom2, atom3, atom4)

    return (length, A1, A2, A3, A4, dihed_angle, x2, y2, z2, x3, y3, z3,
            v21x, v21y, v21z, v23x, v23y, v23z, v34x, v34y, v34z, v32x, v32y, v32z,
            n1x, n1y, n1z, n2x, n2y, n2z)

def dihedral_pot_variables(molecule, atom1, atom2, atom3, atom4, dtype):
    length, A1, A2, A3, A4, dihed_angle = _dihedral_pot_variables(molecule, atom1, atom2, atom3, atom4, dtype)[0:6]
    return length, A1, A2, A3, A4, dihed_angle

def dihedral_force_variables(molecule, atom1, atom2, atom3, atom4, dtype):
    (length, A1, A2, A3, A4, dihed_angle, x2, y2, z2, x3, y3, z3,
     v21x, v21y, v21z, v23x, v23y, v23z, v34x, v34y, v34z, v32x, v32y, v32z,
     n1x, n1y, n1z, n2x, n2y, n2z) = _dihedral_pot_variables(
        molecule, atom1, atom2, atom3, atom4, dtype)
    norm1 = (v21x * v21x + v21y * v21y + v21z * v21z) ** 0.5  # norm of vector v21 (or v12)
    norm2 = (v23x * v23x + v23y * v23y + v23z * v23z) ** 0.5  # norm of vector v23 (or v32)
    norm3 = (v34x * v34x + v34y * v34y + v34z * v34z) ** 0.5  # norm of vector v34 (or v43)

    force = 0
    for i in range(length):
        force += A2[i] * math.sin((A4[i] * dihed_angle) - A3[i]) * A4[i]
    # first derivative of the dihedral potential Vd

    if norm1 == 0 or norm2 == 0:
        f1x = 0
        f1y = 0
        f1z = 0
    else:
        theta1 = math.acos((v21x * v23x + v21y * v23y + v21z * v23z) / (norm1 * norm2))
        M1 = (n1x * n1x + n1y * n1y + n1z * n1z) ** 0.5
        p1x = n1x / M1
        p1y = n1y / M1
        p1z = n1z / M1
        # p1 is the normalization of n1
        u1 = force / (norm1 * math.sin(theta1))
        # force multiplied by the partial derivative of theta1 according to position (x1, y1, z1)
        f1x = p1x * u1
        f1y = p1y * u1
        f1z = p1z * u1

    if norm2 == 0 or norm3 == 0:
        f4x = 0
        f4y = 0
        f4z = 0
    else:
        theta2 = math.acos((v32x * v34x + v32y * v34y + v32z * v34z) / (norm2 * norm3))
        M2 = (n2x * n2x + n2y * n2y + n2z * n2z) ** 0.5
        p2x = n2x / M2
        p2y = n2y / M2
        p2z = n2z / M2
        # p2 is the normalization of n2
        u4 = force / (norm3 * math.sin(theta2))
        # force multiplied by the partial derivative of theta2 according to position (x4, y4, z4)
        f4x = p2x * u4
        f4y = p2y * u4
        f4z = p2z * u4

    ox = (x2 + x3) / 2
    oy = (y2 + y3) / 2
    oz = (z2 + z3) / 2
    # The point o is the center of bond v23
    vo3x = x3 - ox
    vo3y = y3 - oy
    vo3z = z3 - oz
    # vector parting from o to (x3, y3, z3)

    norm_vo3 = vo3x * vo3x + vo3y * vo3y + vo3z * vo3z
    if norm_vo3 == 0:
        f3x = 0
        f3y = 0
        f3z = 0
    else:
        const = -1 / norm_vo3  # inverse square of the norm of vector vo3
        n_vo3_f4x, n_vo3_f4y, n_vo3_f4z = cross_product_3d(molecule, vo3x, f4x, vo3y, f4y, vo3z, f4z)
        n_v34_f4x, n_v34_f4y, n_v34_f4z = cross_product_3d(molecule, v34x, f4x, v34y, f4y, v34z, f4z)
        n_v21_f4x, n_v21_f4y, n_v21_f4z = cross_product_3d(molecule, v21x, f4x, v21y, f4y, v21z, f4z)
        t3x = const * (n_vo3_f4x + 0.5 * n_v34_f4x + 0.5 * n_v21_f4x)
        t3y = const * (n_vo3_f4y - 0.5 * n_v34_f4y - 0.5 * n_v21_f4y)
        t3z = const * (n_vo3_f4z - 0.5 * n_v34_f4z - 0.5 * n_v21_f4z)
        # t3 is the cross product (vo3 by f3)
        f3x, f3y, f3z = cross_product_3d(molecule, t3x, vo3x, t3y, vo3y, t3z, vo3z)

    f2x = -f1x - f3x - f4x
    f2y = -f1y - f3y - f4y
    f2z = -f1z - f3z - f4z
    return f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z, f4x, f4y, f4z

def _angle_dihedral_potential_gaff(molecule, atom1, atom2, atom3, atom4):
    Vd = 0
    dtype = get_dihed_type(molecule, atom1, atom2, atom3, atom4)
    length, A1, A2, A3, A4, dihed_angle = dihedral_pot_variables(
        molecule, atom1, atom2, atom3, atom4, dtype)
    for i in range(length):
        Vd += A2[i] * (1 + math.cos(A4[i] * dihed_angle - A3[i]))
    return Vd, dihed_angle

def angle_dihedral_potential_gaff(molecule, atom1, atom2, atom3, atom4):
    return _angle_dihedral_potential_gaff(molecule, atom1, atom2, atom3, atom4)[0]

def _angle_dihedral_potential_opls(molecule, atom1, atom2, atom3, atom4):
    Vd = 0
    dtype = get_dihed_type(molecule, atom1, atom2, atom3, atom4)
    length, A1, A2, A3, A4, dihed_angle = dihedral_pot_variables(
        molecule, atom1, atom2, atom3, atom4, dtype)
    for i in range(length):
        Vd += 0.5 * (
                A1[i] * (1 + math.cos(dihed_angle)) + A2[i] * (1 - math.cos(2 * dihed_angle)) +
                A3[i] * (1 + math.cos(3 * dihed_angle)) + A4[i])
    return Vd, dihed_angle

def angle_dihedral_potential_opls(molecule, atom1, atom2, atom3, atom4):
    return _angle_dihedral_potential_opls(molecule, atom1, atom2, atom3, atom4)[0]

def _dihedral_potential(molecule, ForceField='gaff'):
    ndihed = molecule.topology.num_dihedrals
    Vd = 0
    dihed_angle = 0
    if ForceField == 'opls':
        for i in range(0, ndihed*4, 4):
            atom1 = molecule.topology.dihedral_list[i] - 1
            atom2 = molecule.topology.dihedral_list[i + 1] - 1
            atom3 = molecule.topology.dihedral_list[i + 2] - 1
            atom4 = molecule.topology.dihedral_list[i + 3] - 1
            if molecule.topology.type[atom1] == 'hc' or molecule.topology.type[atom4] == 'hc':
                continue
            aux1, aux2 = _angle_dihedral_potential_opls(molecule, atom1, atom2, atom3, atom4)
            Vd += aux1
            dihed_angle = aux2
    else:
        for i in range(0, ndihed*4, 4):
            atom1 = molecule.topology.dihedral_list[i] - 1
            atom2 = molecule.topology.dihedral_list[i + 1] - 1
            atom3 = molecule.topology.dihedral_list[i + 2] - 1
            atom4 = molecule.topology.dihedral_list[i + 3] - 1
            if molecule.topology.type[atom1] == 'hc' or molecule.topology.type[atom4] == 'hc':
                continue
            aux1, aux2 = _angle_dihedral_potential_gaff(molecule, atom1, atom2, atom3, atom4)
            Vd += aux1
            dihed_angle = aux2

    return Vd, dihed_angle

def dihedral_potential(molecule, ForceField='gaff'):
    return _dihedral_potential(molecule, ForceField)[0]

def dihedral_force(molecule, fx, fy, fz):
    for i in range(molecule.topology.num_dihedrals):
        atom1 = molecule.topology.dihedral_list[i] - 1
        atom2 = molecule.topology.dihedral_list[i + 1] - 1
        atom3 = molecule.topology.dihedral_list[i + 2] - 1
        atom4 = molecule.topology.dihedral_list[i + 3] - 1
        dtype = get_dihed_type(molecule, atom1, atom2, atom3, atom4)
        if 'hc' in dtype:
            continue
        f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z, f4x, f4y, f4z = dihedral_force_variables(
            molecule, atom1, atom2, atom3, atom4, dtype)
        fx[atom1] += f1x
        fy[atom1] += f1y
        fz[atom1] += f1z
        fx[atom2] += f2x
        fy[atom2] += f2y
        fz[atom2] += f2z
        fx[atom3] += f3x
        fy[atom3] += f3y
        fz[atom3] += f3z
        fx[atom4] += f4x
        fy[atom4] += f4y
        fz[atom4] += f4z
# %%

##############################################################################
if is_main():

    molecule = Molecule()

    path = './butane'
    path1 = path + '/sqm/sqm'
    path2 = path + '/butane'
    '''
    
    path = './3-metil-pentano'
    path1 = path + '/3-metil-pentano'
    path2 = path + '/3-metil-pentano'

    '''
    molecule.read_mol2(path1 + '.mol2')
    molecule.gen_dihed_list_from_angle_list()
    molecule.read_frcmod(path2 + '.frcmod')

    Vb = bond_potential(molecule)
    Va = angle_potential(molecule)
    Vd = dihedral_potential(molecule)

    # minimize_pot_energy(molecule, fx, fy, fz)

    # heuristic_conformation(molecule, theta=np.pi/6, arquivo=path1+'_heuristic.mol2')

    Ud, dihed_angles = dihedral_Ep_rotate(molecule, 2, 3, np.pi/180, 360, 'Ep_sqm2.pdf', False)
    Vb1 = bond_potential(molecule)
    Va1 = angle_potential(molecule)
    Vd1 = dihedral_potential(molecule)
