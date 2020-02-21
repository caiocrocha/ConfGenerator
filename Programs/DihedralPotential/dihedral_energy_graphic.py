import numpy as np
import matplotlib.pyplot as plt

import rotate
import DihedralPotential.dihedral_potential_energy as dihedral_potential_energy
import DihedralPotential.plot_PE as plot_PE
import DihedralPotential.get_dihedral_angle as get_dihedral_angle
import get_bond_vector

def dihedral_energy_graphic(molecule, a, b, c, d, theta, ntimes, pdf_name,
                       write_mol2=False, mol2_name=None):
    times = ntimes + 1
    Vd = np.zeros(times, dtype='float32')
    dihed_angle = np.zeros(times, dtype='float32')
    Vd[0] = dihedral_potential_energy.total_dihedral_potential(molecule)
    v21x, v21y, v21z = get_bond_vector.get_bond_vector_v12(molecule, b-1, a-1)
    v23x, v23y, v23z = get_bond_vector.get_bond_vector_v12(molecule, b-1, c-1)
    v34x, v34y, v34z = get_bond_vector.get_bond_vector_v12(molecule, c-1, d-1)
    angle_rad = get_dihedral_angle.get_dihedral_angle(v21x, v21y, v21z, v23x, v23y, v23z, v34x, v34y, v34z)
    dihed_angle[0] = angle_rad * 180 / np.pi
    plt.scatter(dihed_angle[0], Vd[0], marker='.', color='royalblue')
    rotate.init_rotation_axis(molecule, b-1, c-1, theta)
    rlist = rotate.get_rotation_list(molecule, b-1, c-1)
    if not write_mol2:
        for i in range(1, times):
            for atom in rlist:
                rotate.rotate_atom(molecule, atom)
            Vd[i] = dihedral_potential_energy.total_dihedral_potential(molecule)
            v21x, v21y, v21z = get_bond_vector.get_bond_vector_v12(molecule, b-1, a-1)
            v23x, v23y, v23z = get_bond_vector.get_bond_vector_v12(molecule, b-1, c-1)
            v34x, v34y, v34z = get_bond_vector.get_bond_vector_v12(molecule, c-1, d-1)
            angle_rad = get_dihedral_angle.get_dihedral_angle(v21x, v21y, v21z, v23x, v23y, v23z, v34x, v34y, v34z)
            dihed_angle[i] = angle_rad * 180 / np.pi
            plt.scatter(dihed_angle[i], Vd[i], marker='.', color='royalblue')
    else:
        molecule.write_mol2(mol2_name, 'w+')
        for i in range(1, times):
            for atom in rlist:
                rotate.rotate_atom(molecule, atom)
            Vd[i] = dihedral_potential_energy.total_dihedral_potential(molecule)
            v21x, v21y, v21z = get_bond_vector.get_bond_vector_v12(molecule, b-1, a-1)
            v23x, v23y, v23z = get_bond_vector.get_bond_vector_v12(molecule, b-1, c-1)
            v34x, v34y, v34z = get_bond_vector.get_bond_vector_v12(molecule, c-1, d-1)
            angle_rad = get_dihedral_angle.get_dihedral_angle(v21x, v21y, v21z, v23x, v23y, v23z, v34x, v34y, v34z)
            dihed_angle[i] = angle_rad * 180 / np.pi
            plt.scatter(dihed_angle[i], Vd[i], marker='.', color='royalblue')
            molecule.write_mol2(mol2_name, 'a')
    plot_PE.plot_PE(theta, ntimes, pdf_name)
    return dihed_angle, Vd