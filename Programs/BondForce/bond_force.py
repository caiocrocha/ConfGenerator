import BondPotential.bond_potential_energy as bond_potential_energy

def bond_force_variables(molecule, atom1, atom2):
    bx, by, bz, b, b0, Kb, d = bond_potential_energy.bond_pot_variables(molecule, atom1, atom2)
    if b != 0:
        ux = bx / b
        uy = by / b
        uz = bz / b
    else:
        ux = 0
        uy = 0
        uz = 0
    return bx, by, bz, b, b0, Kb, d, ux, uy, uz

def bond_force(molecule, fx, fy, fz):
    for i in range(0, molecule.topology.num_bonds * 2, 2):
        atom1 = molecule.topology.bond_list[i] - 1
        atom2 = molecule.topology.bond_list[i + 1] - 1
        bx, by, bz, b, b0, Kb, d, ux, uy, uz = bond_force_variables(molecule, atom1, atom2)
        force = -2 * Kb * d  # first derivative of the bond potential Vb
        fx1 = force * ux  # decomposition of the force on the X axis
        fy1 = force * uy  # decomposition of the force on the Y axis
        fz1 = force * uz  # decomposition of the force on the Z axis
        fx2 = -fx1
        fy2 = -fy1
        fz2 = -fz1
        fx[atom1] += fx1
        fy[atom1] += fy1
        fz[atom1] += fz1
        fx[atom2] += fx2
        fy[atom2] += fy2
        fz[atom2] += fz2