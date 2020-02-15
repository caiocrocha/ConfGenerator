def cross_product_3d(v1x, v2x, v1y, v2y, v1z, v2z):
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