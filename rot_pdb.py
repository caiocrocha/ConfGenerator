#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 10:15:23 2019

@author: dgomes
"""

def rot_pdb(coordenadas, a, b, theta):
    
    
    x=coordenadas.x
    y=coordenadas.y
    z=coordenadas.z
    
    
    # Calcula quais os atomos que eu preciso rodar
    # A ordem do a,b importa !!!
    
    
    for i in range(len(coordenadas)):
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
return(coordenadas)