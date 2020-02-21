def get_torsions_list(molecule, atom1, atom2, atom3, atom4):
    dihedrals_list = []
    k = int(atom3 * (atom3 - 1) / 2)

    for i in range(atom3):
        if (molecule.topology.bond_matrix[k + i] and i != atom2):
            dihedrals_list.append(atom1 + 1)
            dihedrals_list.append(atom2 + 1)
            dihedrals_list.append(atom3 + 1)
            dihedrals_list.append(i + 1)

    for i in range(atom3 + 1, molecule.num_atoms):
        if (molecule.topology.bond_matrix[int(i * (i - 1) / 2 + atom3)] and i != atom2):
            dihedrals_list.append(atom1 + 1)
            dihedrals_list.append(atom2 + 1)
            dihedrals_list.append(atom3 + 1)
            dihedrals_list.append(i + 1)

    return dihedrals_list
