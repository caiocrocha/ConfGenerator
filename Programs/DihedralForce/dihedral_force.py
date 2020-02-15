import math

import cross_product_3d
import DihedralPotential.dihedral_potential_energy as dihedral_potential_energy
import get_dihedral_type

def dihedral_force_variables(molecule, atom1, atom2, atom3, atom4, dtype):
    (length, A1, A2, A3, A4, dihed_angle, x2, y2, z2, x3, y3, z3,
     v21x, v21y, v21z, v23x, v23y, v23z, v34x, v34y, v34z, v32x, v32y, v32z,
     n1x, n1y, n1z, n2x, n2y, n2z) = dihedral_potential_energy._dihedral_pot_variables(
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
        n_vo3_f4x, n_vo3_f4y, n_vo3_f4z = cross_product_3d.cross_product_3d(molecule, vo3x, f4x, vo3y, f4y, vo3z, f4z)
        n_v34_f4x, n_v34_f4y, n_v34_f4z = cross_product_3d.cross_product_3d(molecule, v34x, f4x, v34y, f4y, v34z, f4z)
        n_v21_f4x, n_v21_f4y, n_v21_f4z = cross_product_3d.cross_product_3d(molecule, v21x, f4x, v21y, f4y, v21z, f4z)
        t3x = const * (n_vo3_f4x + 0.5 * n_v34_f4x + 0.5 * n_v21_f4x)
        t3y = const * (n_vo3_f4y - 0.5 * n_v34_f4y - 0.5 * n_v21_f4y)
        t3z = const * (n_vo3_f4z - 0.5 * n_v34_f4z - 0.5 * n_v21_f4z)
        # t3 is the cross product (vo3 by f3)
        f3x, f3y, f3z = cross_product_3d.cross_product_3d(molecule, t3x, vo3x, t3y, vo3y, t3z, vo3z)

    f2x = -f1x - f3x - f4x
    f2y = -f1y - f3y - f4y
    f2z = -f1z - f3z - f4z
    return f1x, f1y, f1z, f2x, f2y, f2z, f3x, f3y, f3z, f4x, f4y, f4z

def dihedral_force(molecule, fx, fy, fz):
    for i in range(molecule.topology.num_dihedrals):
        atom1 = molecule.topology.dihedral_list[i] - 1
        atom2 = molecule.topology.dihedral_list[i + 1] - 1
        atom3 = molecule.topology.dihedral_list[i + 2] - 1
        atom4 = molecule.topology.dihedral_list[i + 3] - 1
        dtype = get_dihedral_type.get_dihedral_type(molecule, atom1, atom2, atom3, atom4)
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