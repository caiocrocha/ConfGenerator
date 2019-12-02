#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 09:51:57 2019

@author: dgomes
"""

import numpy as np
import rotate

natm=4
a=0
b=1
r=[2, 3]
theta = np.pi/2

x = [0.0, 0.0, 1.0, 1.0]
y = [0.0, -1.0, 0.0, 1.0]
z = [0.0, 0.0, 0.0, 0.0]

for i in range(natm):
    print(x[i], y[i], z[i])

Ax = x[a]
Ay = y[a]
Az = z[a]
Bx = x[b]
By = y[b]
Bz = z[b]
M = (((Bx - Ax)**2) + ((By - Ay)**2) + ((Bz - Az)**2))**0.5
# M is the magnitude of the unit vector u
ux = (Bx - Ax)/M # x component of u
uy = (By - Ay)/M # y component of u
uz = (Bz - Az)/M # z component of u
cos0 = np.cos(theta) # cosine of the angle of rotation
sin0 = np.sin(theta) # sine of the angle of rotation
t = 1 - cos0
for i in range(len(r)):
    j = r[i]
    x1 = x[j] - Ax
    y1 = y[j] - Ay
    z1 = z[j] - Az
    x[j] = round(((cos0 + t*(ux**2))*x1 + (t*ux*uy - sin0*uz)*y1 + (t*ux*uz + sin0*uy)*z1), 3) + Ax
    y[j] = round(((t*ux*uy + sin0*uz)*x1 + (cos0 + t*(uy**2))*y1 + (t*uy*uz - sin0*ux)*z1), 3) + Ay
    z[j] = round(((t*ux*uz - sin0*uy)*x1 + (t*uy*uz + sin0*ux)*y1 + (cos0 + t*(uz**2))*z1), 3) + Az
    # rotation matrix applied to point
    
print("Matriz depois da rotacao: ")

for i in range(4):
    print(x[i], y[i], z[i])