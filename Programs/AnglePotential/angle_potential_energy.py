import math

def _angle_pot_variables(molecule, atom1, atom2, atom3):

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

    atype = ('{}-{}-{}'.format(molecule.atom_type[atom1],
                               molecule.atom_type[atom2], molecule.atom_type[atom3])
             )  # angle type (e.g. 'c3-c3-hc' or 'c3-c3-c3')
    try:
        Ka = molecule.topology.angle_types[atype][0]  # angle elastic constant
        a0 = molecule.topology.angle_types[atype][1]  # equilibrium angle
    except:
        atype = ('{}-{}-{}'.format(molecule.atom_type[atom3],
                                   molecule.atom_type[atom2], molecule.atom_type[atom1])
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