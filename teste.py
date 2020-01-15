import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

#%%
def ang2grad(a) :
    return a*np.pi/180

A2=[0.18,0.25,0.20]
A3=[0,180,180]
A4=[-3,-2,1]

data=[]
for ang in np.linspace(0,360,361) :
    Vd=0
    for j in range(len(A2)) : 
        Vd += A2[j] * (1 + math.cos( A4[j]* ang2grad(ang) - ang2grad(A3[j] ))) 
    
    data.append([ang,Vd])
    
data=pd.DataFrame(data)
data.columns=['ang','Vd']
data.plot(x='ang',y='Vd')
data.to_csv('curva.dat',index=False, float_format='%.3f')
