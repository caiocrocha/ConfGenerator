import math

import rotate
import DihedralPotential.dihedral_angle_potential as dihedral_angle_potential
import DihedralPotential.dihedral_potential_energy as dihedral_potential_energy

def heuristic_rotate(molecule, rlist, theta, ntimes, atom1, atom2, atom3, atom4):
    angles_list = []
    rotate.init_rotation_axis(molecule, atom2, atom3, theta)
    lowest, dihed0 = dihedral_angle_potential._dihedral_angle_potential_gaff(molecule, atom1, atom2, atom3, atom4)
    dihed0 *= 180 / math.pi
    angles_list.append(dihed0)
    aux1 = dihedral_potential_energy.total_dihedral_potential(molecule)
    comeback = 0
    for i in range(ntimes - 1):
        for atom in rlist:
            rotate.rotate_atom(molecule, atom)
        Vd, dihed1 = dihedral_angle_potential._dihedral_angle_potential_gaff(molecule, atom1, atom2, atom3, atom4)
        dihed1 *= 180/math.pi
        angles_list.append(dihed1)
        aux2 = dihedral_potential_energy.total_dihedral_potential(molecule)
        if aux2 < aux1:
            lowest = Vd
            comeback = i + 1
            aux1 = aux2
    rotate.set_rotation_angle(theta * (1 + comeback))
    for atom in rlist:
        rotate.rotate_atom(molecule, atom)

def heuristic_conformation(molecule, theta):
    ntimes = int(2 * math.pi / theta)
    for i in range(0, molecule.topology.num_dihedrals*4, 4):
        atom1 = molecule.topology.dihedral_list[i] - 1
        atom2 = molecule.topology.dihedral_list[i + 1] - 1
        atom3 = molecule.topology.dihedral_list[i + 2] - 1
        atom4 = molecule.topology.dihedral_list[i + 3] - 1
        if molecule.topology.type[atom1] == 'hc' or molecule.topology.type[atom4] == 'hc':
            continue
        rlist = rotate.get_rotation_list(molecule, atom2, atom3)
        heuristic_rotate(molecule, rlist, theta, ntimes, atom1, atom2, atom3, atom4)