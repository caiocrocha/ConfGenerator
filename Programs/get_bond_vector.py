def get_bond_vector_v12(molecule, atom1, atom2):
    return (molecule.x[atom2] - molecule.x[atom1],
            molecule.y[atom2] - molecule.y[atom1],
            molecule.z[atom2] - molecule.z[atom1])