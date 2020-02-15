import numpy as np

import cross_product_3d
import get_dihedral_angle

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

    n1x, n1y, n1z = cross_product_3d.cross_product_3d(v21x, v23x, v21y, v23y, v21z, v23z)
    # normal vector of the plane determined by vectors v21 and v23
    n2x, n2y, n2z = cross_product_3d.cross_product_3d(v32x, v34x, v32y, v34y, v32z, v34z)
    # normal vector of the plane determined by vectors v32 and v34
    '''
    nx, ny, nz = cross_product_3d.cross_product_3d(molecule, n1x, n2x, n1y, n2y, n1z, n2z)
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
    dihed_angle = get_dihedral_angle.get_dihedral_angle(molecule, atom1, atom2, atom3, atom4)

    return (length, A1, A2, A3, A4, dihed_angle, x2, y2, z2, x3, y3, z3,
            v21x, v21y, v21z, v23x, v23y, v23z, v34x, v34y, v34z, v32x, v32y, v32z,
            n1x, n1y, n1z, n2x, n2y, n2z)

def dihedral_pot_variables(molecule, atom1, atom2, atom3, atom4, dtype):
    length, A1, A2, A3, A4, dihed_angle = _dihedral_pot_variables(molecule, atom1, atom2, atom3, atom4, dtype)[0:6]
    return length, A1, A2, A3, A4, dihed_angle