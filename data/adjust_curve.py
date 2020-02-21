import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt

datay = [4.5, 4.0, 3.3, 2.6, 2.1, 1.2, 0.9, 1.2, 1.7, 2.3, 2.8, 3.3, 3.5, 3.3, 2.6, 1.7, 1.1, 0.4, 0.0, 0.4, 1.1, 1.7, 2.6, 3.3, 3.5, 3.3, 2.8, 2.3, 1.7, 1.2, 0.9, 1.2, 2.1, 2.6, 3.3, 4.0, 4.5]
datax = np.linspace(0, 360, 37)

#		            Vd	  lambda	n
#c3-c3-c3-c3	1	0,18	0	  -3
#c3-c3-c3-c3	1	0,25	180	  -2
#c3-c3-c3-c3	1	0,2	    180	   1


def func(ang, vd1,l1,vd2,l2,vd3,l3):
    ang = np.radians(ang)
#    return ( vd1 * ( 1 + np.cos(n1 * ang - l1)) + vd2 * ( 1 + np.cos(n2 * ang - l2) + vd3 * ( 1 + np.cos(n3 * ang - l3) )))
    return ( vd1 * ( 1 + np.cos(-3 * ang - l1)) + vd2 * ( 1 + np.cos(-2 * ang - l2)) + vd3 * ( 1 + np.cos(1 * ang - l3)) )

popt, pcov = scipy.optimize.curve_fit(func, datax, datay)

plt.plot(datax, datay, 'b-', label='data')
plt.plot(datax, func(datax, *popt), 'r-', label='fit: vd1=%5.3f, l1=%5.3f,vd2=%5.3f, l2=%5.3f, vd3=%5.3f, l3=%5.3f' % tuple(popt))
plt.xlabel('Angle')
plt.ylabel('Potential')
plt.legend()
plt.savefig('compare_curves_deg.pdf')
plt.show()

