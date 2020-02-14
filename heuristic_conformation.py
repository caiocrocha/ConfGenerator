import numpy as np

import rotate
import angle_dihedral_potential
import dihedral_potential_energy

def heuristic_rotate(molecule, rlist, theta, ntimes, arquivo, atom1, atom2, atom3, atom4, vez):
    rotate.init_rotation_axis(molecule, atom2 + 1, atom3 + 1, theta)
    lowest, dihed0 = angle_dihedral_potential._angle_dihedral_potential_gaff(molecule, atom1, atom2, atom3, atom4)
    dihed0 *= 180 / np.pi
    aux1 = dihedral_potential_energy.dihedral_potential(molecule)
    comeback = 0
    for i in range(ntimes - 1):
        for atom in rlist:
            rotate.rotate_atom(molecule, atom)
        # molecule.write_pdb(arquivo, 'a', vez)
        vez += 1
        Vd, dihed1 = angle_dihedral_potential._angle_dihedral_potential_gaff(molecule, atom1, atom2, atom3, atom4)
        dihed1 *= 180/np.pi
        aux2 = dihedral_potential_energy.dihedral_potential(molecule)
        if aux2 < aux1:
            lowest = Vd
            comeback = i + 1
            aux1 = aux2
    rotate.set_rotation_angle(theta * (1 + comeback))
    for atom in rlist:
        rotate.rotate_atom(molecule, atom)
    # molecule.write_pdb(arquivo, 'a', vez)
    vez += 1

def heuristic_conformation(molecule, theta, arquivo):
    ntimes = int(2 * np.pi / theta)
    global vez
    vez = 1
    # molecule.write_pdb(arquivo, 'w+', 0)
    for i in range(0, molecule.topology.num_dihedrals, 4):
        atom1 = molecule.topology.dihedral_list[i] - 1
        atom2 = molecule.topology.dihedral_list[i + 1] - 1
        atom3 = molecule.topology.dihedral_list[i + 2] - 1
        atom4 = molecule.topology.dihedral_list[i + 3] - 1
        if molecule.topology.type[atom1] == 'hc' or molecule.topology.type[atom4] == 'hc':
            continue
        rlist = rotate.get_rotation_list(molecule, atom2 + 1, atom3 + 1)
        heuristic_rotate(molecule, rlist, theta, ntimes, arquivo, atom1, atom2, atom3, atom4, vez)
    molecule.write_mol2(arquivo)