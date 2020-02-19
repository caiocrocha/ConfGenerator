import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt

datay = [4.5, 4.0, 3.3, 2.6, 2.1, 1.2, 0.9, 1.2, 1.7, 2.3, 2.8, 3.3, 3.5, 3.3, 2.6, 1.7, 1.1, 0.4, 0.0, 0.4, 1.1, 1.7, 2.6, 3.3, 3.5, 3.3, 2.8, 2.3, 1.7, 1.2, 0.9, 1.2, 2.1, 2.6, 3.3, 4.0, 4.5]
datax = np.linspace(0, np.pi*2, 37)

def func(x, a, b, c):
    return a * (1 + np.cos(c * x - b))

popt, pcov = scipy.optimize.curve_fit(func, datax, datay)

plt.plot(datax, datay, 'b-', label='data')
plt.plot(datax, func(datax, *popt), 'r-', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
plt.xlabel('Angle')
plt.ylabel('Potential')
plt.legend()
plt.savefig('compare_curves_rad.pdf')
plt.show()
