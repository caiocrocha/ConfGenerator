import numpy as np
import matplotlib.pyplot as plt

import rotate

def dihedral_energy_graphic(molecule, a, b, theta, ntimes, pdf_name,
                       write_pdb=False, pdb_name=None):
    times = ntimes + 1
    Vd = np.zeros(times, dtype='float32')
    dihed_angle = np.zeros(times, dtype='float32')
    Vd[0], dihed_angle[0] = dihedral_potential_energy._dihedral_potential(molecule)
    plt.scatter(dihed_angle[0], Vd[0], marker='.', color='royalblue')
    rotate.init_rotation_axis(molecule, a, b, theta)
    rlist = rotate.get_rotation_list(molecule, a, b)
    if not write_pdb:
        for i in range(1, times):
            for atom in rlist:
                rotate.rotate_atom(molecule, atom)
            Vd[i], dihed_angle[i] = dihedral_potential_energy._dihedral_potential(molecule)
            dihed_angle[i] = dihed_angle[i]*180/np.pi + 180
            plt.scatter(dihed_angle[i], Vd[i], marker='.', color='royalblue')
    else:
        molecule.write_pdb(pdb_name, 'w+', 0)
        for i in range(1, times):
            for atom in rlist:
                rotate.rotate_atom(molecule, atom)
            Vd[i], dihed_angle[i] = dihedral_potential_energy._dihedral_potential(molecule)
            dihed_angle[i] = dihed_angle[i]*180/np.pi + 180
            plt.scatter(dihed_angle[i], Vd[i], marker='.', color='royalblue')
            molecule.write_pdb(pdb_name, 'a', i)
    plot_PE.plot_PE(theta, ntimes, pdf_name)
    return dihed_angle, Vd