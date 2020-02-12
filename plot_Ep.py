import math
import matplotlib.pyplot as plt

def plot_Ep(theta, ntimes, pdf_name):
    plt.xlabel('Degrees (rad)')
    plt.ylabel(r'Elastic potential (kcal $mol^{-1} \AA^{-2}$)')
    plt.grid()
    plt.suptitle('Rotation degrees X Dihedral elastic potential (Vd) graphic')
    plt.title(r'{}$^\circ$ x {} rotations'.format(theta * 180 / math.pi, ntimes), fontsize=10, loc='right')
    plt.savefig(pdf_name)
    plt.show()