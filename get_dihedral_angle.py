import numpy as np

def get_dihedral_angle(molecule, atom1, atom2, atom3, atom4):
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