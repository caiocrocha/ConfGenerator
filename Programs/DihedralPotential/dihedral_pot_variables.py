import numpy as np

import cross_product_3d
import DihedralPotential.get_dihedral_angle as get_dihedral_angle
import DihedralPotential.get_dihedral_type as get_dihedral_type

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

    dihed_angle = get_dihedral_angle.get_dihedral_angle(v21x, v21y, v21z, v23x, v23y, v23z, v34x, v34y, v34z)

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
        dtype = get_dihedral_type.get_dihedral_type(molecule, atom4, atom3, atom2, atom1)
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

    return (length, A1, A2, A3, A4, dihed_angle, x2, y2, z2, x3, y3, z3,
            v21x, v21y, v21z, v23x, v23y, v23z, v34x, v34y, v34z, v32x, v32y, v32z)

def dihedral_pot_variables(molecule, atom1, atom2, atom3, atom4, dtype):
    length, A1, A2, A3, A4, dihed_angle = _dihedral_pot_variables(molecule, atom1, atom2, atom3, atom4, dtype)[0:6]
    return length, A1, A2, A3, A4, dihed_angle