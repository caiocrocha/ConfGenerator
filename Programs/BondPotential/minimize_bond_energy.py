def minimize_bond_energy(molecule, fx, fy, fz):
    for i in range(molecule.num_atoms):
        fx1 = fx[i]
        fy1 = fy[i]
        fz1 = fz[i]
        N = (fx1 * fx1 + fy1 * fy1 + fz1 * fz1) ** 0.5  # normal force
        molecule.x[i] -= 0.01 * fx1 / N
        molecule.y[i] -= 0.01 * fy1 / N
        molecule.z[i] -= 0.01 * fz1 / N