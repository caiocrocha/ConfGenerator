import numpy as np
import matplotlib.pyplot as plt

import rotate
import DihedralPotential.dihedral_potential_energy as dihedral_potential_energy
import DihedralPotential.plot_PE as plot_PE

def dihedral_energy_graphic(molecule, a, b, theta, ntimes, pdf_name,
                       write_mol2=False, mol2_name=None):
    times = ntimes + 1
    Vd = np.zeros(times, dtype='float32')
    dihed_angle = np.zeros(times, dtype='float32')
    Vd[0], dihed_angle[0] = dihedral_potential_energy._total_dihedral_potential(molecule)
    plt.scatter(dihed_angle[0], Vd[0], marker='.', color='royalblue')
    rotate.init_rotation_axis(molecule, a-1, b-1, theta)
    rlist = rotate.get_rotation_list(molecule, a-1, b-1)
    if not write_mol2:
        for i in range(1, times):
            for atom in rlist:
                rotate.rotate_atom(molecule, atom)
            Vd[i], dihed_angle[i] = dihedral_potential_energy._total_dihedral_potential(molecule)
            dihed_angle[i] = dihed_angle[i]*180/np.pi + 180
            plt.scatter(dihed_angle[i], Vd[i], marker='.', color='royalblue')
    else:
        molecule.write_mol2(mol2_name + '0.mol2')
        for i in range(1, times):
            for atom in rlist:
                rotate.rotate_atom(molecule, atom)
            Vd[i], dihed_angle[i] = dihedral_potential_energy._total_dihedral_potential(molecule)
            dihed_angle[i] = dihed_angle[i]*180/np.pi + 180
            plt.scatter(dihed_angle[i], Vd[i], marker='.', color='royalblue')
            molecule.write_mol2(mol2_name + '{}.mol2'.format(i))
    plot_PE.plot_PE(theta, ntimes, pdf_name)
    return dihed_angle, Vd