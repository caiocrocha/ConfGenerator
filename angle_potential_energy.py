import angle_pot_variables

def angle_potential(molecule):
    nangles = molecule.topology.num_angles
    Va = 0
    for i in range(0, nangles * 3, 3):
        delta_angle, Ka = angle_pot_variables.angle_pot_variables(
            molecule,
            atom1=molecule.topology.angle_list[i] - 1,
            atom2=molecule.topology.angle_list[i + 1] - 1,
            atom3=molecule.topology.angle_list[i + 2] - 1)
        Va += Ka * delta_angle * delta_angle
    return Va