def get_dihedral_type(molecule, atom1, atom2, atom3, atom4):
    dtype = ('{}-{}-{}-{}'.format(molecule.atom_type[atom1],
                                  molecule.atom_type[atom2],
                                  molecule.atom_type[atom3],
                                  molecule.atom_type[atom4])
             )  # dihedral type (e.g. 'c3-c3-c3-c3')
    return dtype