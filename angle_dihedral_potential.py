import math

import get_dihedral_type
import dihedral_pot_variables

def _angle_dihedral_potential_gaff(molecule, atom1, atom2, atom3, atom4):
    Vd = 0
    dtype = get_dihedral_type.get_dihedral_type(molecule, atom1, atom2, atom3, atom4)
    length, A1, A2, A3, A4, dihed_angle = dihedral_pot_variables.dihedral_pot_variables(
        molecule, atom1, atom2, atom3, atom4, dtype)
    for i in range(length):
        Vd += A2[i] * (1 + math.cos(A4[i] * dihed_angle - A3[i]))
    return Vd, dihed_angle

def angle_dihedral_potential_gaff(molecule, atom1, atom2, atom3, atom4):
    return _angle_dihedral_potential_gaff(molecule, atom1, atom2, atom3, atom4)[0]

def _angle_dihedral_potential_opls(molecule, atom1, atom2, atom3, atom4):
    Vd = 0
    dtype = get_dihedral_type.get_dihed_type(molecule, atom1, atom2, atom3, atom4)
    length, A1, A2, A3, A4, dihed_angle = dihedral_pot_variables.dihedral_pot_variables(
        molecule, atom1, atom2, atom3, atom4, dtype)
    for i in range(length):
        Vd += 0.5 * (
                A1[i] * (1 + math.cos(dihed_angle)) + A2[i] * (1 - math.cos(2 * dihed_angle)) +
                A3[i] * (1 + math.cos(3 * dihed_angle)) + A4[i])
    return Vd, dihed_angle

def angle_dihedral_potential_opls(molecule, atom1, atom2, atom3, atom4):
    return _angle_dihedral_potential_opls(molecule, atom1, atom2, atom3, atom4)[0]