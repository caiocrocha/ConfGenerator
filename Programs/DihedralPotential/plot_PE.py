import math
import matplotlib.pyplot as plt

def plot_PE(theta, ntimes, pdf_name):
    plt.xlabel('Dihedral angle (degrees)')
    plt.ylabel(r'Elastic potential (kcal $mol^{-1} \AA^{-2}$)')
    plt.grid()
    plt.suptitle('Rotation Degrees x Dihedral Elastic Potential (Vd)')
    plt.title(r'{}$^\circ$ x {} rotations'.format(theta * 180 / math.pi, ntimes), fontsize=10, loc='right')
    plt.savefig(pdf_name)
    plt.show()
