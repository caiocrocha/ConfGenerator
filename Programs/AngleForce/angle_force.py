import math

import cross_product_3d

def angle_force_variables(molecule, atom1, atom2, atom3):
    a0, Ka, v21x, v23x, v21y, v23y, v21z, v23z, norm1, norm2 = angle_pot_variables._angle_pot_variables(
        molecule, atom1, atom2, atom3)
    nx, ny, nz = cross_product_3d.cross_product_3d(v21x, v23x, v21y, v23y, v21z, v23z)
    # normal vector of the plane determined by vectors v21 and v23
    vR1x, vR1y, vR1z = cross_product_3d.cross_product_3d(v21x, nx, v21y, ny, v21z, nz)
    # vR1: resulting vector (orthogonal to v21)
    vR2x, vR2y, vR2z = cross_product_3d.cross_product_3d(v23x, nx, v23y, ny, v23z, nz)
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
