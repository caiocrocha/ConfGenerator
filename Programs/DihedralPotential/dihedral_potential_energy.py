import DihedralPotential.dihedral_angle_potential as dihedral_angle_potential

def _total_dihedral_potential(molecule, ForceField='gaff'):
    ndihed = molecule.topology.num_dihedrals
    Vd = 0
    dihed_angle = 0
    if ForceField == 'opls':
        for i in range(0, ndihed*4, 4):
            atom1 = molecule.topology.dihedral_list[i] - 1
            atom2 = molecule.topology.dihedral_list[i + 1] - 1
            atom3 = molecule.topology.dihedral_list[i + 2] - 1
            atom4 = molecule.topology.dihedral_list[i + 3] - 1
            if molecule.topology.type[atom1] == 'hc' or molecule.topology.type[atom4] == 'hc':
                continue
            aux1, aux2 = dihedral_angle_potential._dihedral_angle_potential_gaff(molecule, atom1, atom2, atom3, atom4)
            Vd += aux1
            dihed_angle = aux2
    else:
        for i in range(0, ndihed*4, 4):
            atom1 = molecule.topology.dihedral_list[i] - 1
            atom2 = molecule.topology.dihedral_list[i + 1] - 1
            atom3 = molecule.topology.dihedral_list[i + 2] - 1
            atom4 = molecule.topology.dihedral_list[i + 3] - 1
            if molecule.topology.type[atom1] == 'hc' or molecule.topology.type[atom4] == 'hc':
                continue
            aux1, aux2 = dihedral_angle_potential._dihedral_angle_potential_gaff(molecule, atom1, atom2, atom3, atom4)
            Vd += aux1
            dihed_angle = aux2

    return Vd, dihed_angle

def total_dihedral_potential(molecule, ForceField='gaff'):
    return _total_dihedral_potential(molecule, ForceField)[0]