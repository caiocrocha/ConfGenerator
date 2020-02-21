def bond_pot_variables(molecule, atom1, atom2):
    bx = molecule.x[atom2] - molecule.x[atom1]
    by = molecule.y[atom2] - molecule.y[atom1]
    bz = molecule.z[atom2] - molecule.z[atom1]
    b = (bx * bx + by * by + bz * bz) ** 0.5  # distance between atoms 1 and 2
    btype = ('{}-{}'.format(molecule.atom_type[atom1],
                            molecule.atom_type[atom2])
             )  # bond type (e.g. 'c3-c3' or 'c3-hc')
    try:
        Kb = molecule.topology.bond_types[btype][0]  # elastic constant of the bond
        b0 = molecule.topology.bond_types[btype][1]  # equilibrium bond length
    except:
        btype = ('{}-{}'.format(molecule.atom_type[atom2],
                                molecule.atom_type[atom1])
                 )
        try:
            Kb = molecule.topology.bond_types[btype][0]
            b0 = molecule.topology.bond_types[btype][1]
        except:
            pass
    d = b - b0
    # d: difference between the present distance between the two atoms
    # and the equilibrium distance
    return bx, by, bz, b, b0, Kb, d

def bond_potential(molecule):
    nbonds = molecule.topology.num_bonds
    Vb = 0
    for i in range(0, nbonds * 2, 2):
        bx, by, bz, b, b0, Kb, d = bond_pot_variables(
            molecule,
            atom1=molecule.topology.bond_list[i] - 1,
            atom2=molecule.topology.bond_list[i + 1] - 1)
        Vb += Kb * d * d
    return Vb
