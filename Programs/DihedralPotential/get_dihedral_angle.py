import math

import cross_product_3d

def get_dihedral_angle(v21x, v21y, v21z, v23x, v23y, v23z, v34x, v34y, v34z):
    # Adapted from https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python

    norm23 = (v23x*v23x + v23y*v23y + v23z*v23z)**0.5
    v23x /= norm23
    v23y /= norm23
    v23z /= norm23

    dot1 = v21x * v23x + v21y * v23y + v21z * v23z
    proj1x = v21x - dot1 * v23x
    proj1y = v21y - dot1 * v23y
    proj1z = v21z - dot1 * v23z

    dot2 = v34x * v23x + v34y * v23y + v34z * v23z
    proj2x = v34x - dot2 * v23x
    proj2y = v34y - dot2 * v23y
    proj2z = v34z - dot2 * v23z

    dx = proj1x*proj2x + proj1y*proj2y + proj1z*proj2z
    cx, cy, cz = cross_product_3d.cross_product_3d(v23x, proj1x, v23y, proj1y, v23z, proj1z)
    dy = cx*proj2x + cy*proj2y + cz*proj2z

    return math.atan2(dy, dx)
