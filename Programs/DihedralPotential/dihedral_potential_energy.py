import DihedralPotential.dihedral_angle_potential as dihedral_angle_potential

def total_dihedral_potential(molecule, ForceField='gaff'):
    ndihed = molecule.topology.num_dihedrals
    Vd = 0
    if ForceField == 'opls':
        for i in range(0, ndihed*4, 4):
            atom1 = molecule.topology.dihedral_list[i] - 1
            atom2 = molecule.topology.dihedral_list[i + 1] - 1
            atom3 = molecule.topology.dihedral_list[i + 2] - 1
            atom4 = molecule.topology.dihedral_list[i + 3] - 1
            if molecule.atom_type[atom1] == 'hc' or molecule.atom_type[atom4] == 'hc':
                continue
            Vd += dihedral_angle_potential.dihedral_angle_potential_gaff(molecule, atom1, atom2, atom3, atom4)
    else:
        for i in range(0, ndihed*4, 4):
            atom1 = molecule.topology.dihedral_list[i] - 1
            atom2 = molecule.topology.dihedral_list[i + 1] - 1
            atom3 = molecule.topology.dihedral_list[i + 2] - 1
            atom4 = molecule.topology.dihedral_list[i + 3] - 1
            if molecule.atom_type[atom1] == 'hc' or molecule.atom_type[atom4] == 'hc':
                continue
            Vd += dihedral_angle_potential.dihedral_angle_potential_gaff(molecule, atom1, atom2, atom3, atom4)

    return Vd

def partial_dihedral_potential(molecule, dihedrals_list, ForceField='gaff'):
    Vd = 0
    if ForceField == 'opls':
        for i in range(0, len(dihedrals_list), 4):
            atom1 = dihedrals_list[i] - 1
            atom2 = dihedrals_list[i + 1] - 1
            atom3 = dihedrals_list[i + 2] - 1
            atom4 = dihedrals_list[i + 3] - 1
            if molecule.atom_type[atom1] == 'hc' or molecule.atom_type[atom4] == 'hc':
                continue
            Vd += dihedral_angle_potential.dihedral_angle_potential_opls(molecule, atom1, atom2, atom3, atom4)
    else:
        for i in range(0, len(dihedrals_list), 4):
            atom1 = dihedrals_list[i] - 1
            atom2 = dihedrals_list[i + 1] - 1
            atom3 = dihedrals_list[i + 2] - 1
            atom4 = dihedrals_list[i + 3] - 1
            if molecule.atom_type[atom1] == 'hc' or molecule.atom_type[atom4] == 'hc':
                continue
            Vd += dihedral_angle_potential.dihedral_angle_potential_gaff(molecule, atom1, atom2, atom3, atom4)

    return Vd
