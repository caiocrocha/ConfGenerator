def get_dihedral_type(molecule, atom1, atom2, atom3, atom4):
    dtype = ('{}-{}-{}-{}'.format(molecule.topology.type[atom1],
                                  molecule.topology.type[atom2],
                                  molecule.topology.type[atom3],
                                  molecule.topology.type[atom4])
             )  # dihedral type (e.g. 'c3-c3-c3-c3')
    return dtype