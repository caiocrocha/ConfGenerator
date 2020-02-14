import numpy as np
import os

from class_molecule import Molecule
import dihedral_energy_graphic
import heuristic_conformation
import angle_potential_energy
import angle_dihedral_potential
import bond_potential_energy
import dihedral_potential_energy
import minimize_bond_energy
import bond_force
import angle_force
import dihedral_force

def is_main():
    return __name__ == '__main__'

##############################################################################
if is_main():

    molecule = Molecule()

    path = os.getcwd()

    path1 = path + '/butane'
    path2 = path1 + '/sqm/sqm'
    path3 = path1 + '/butane'
    '''
    
    path1 = path + '/3-metil-pentano'
    path2 = path1 + '/3-metil-pentano'
    path3 = path1 + '/3-metil-pentano'

    '''
    molecule.read_mol2(path2 + '.mol2')
    molecule.gen_dihed_list_from_angle_list()
    molecule.read_frcmod(path3 + '.frcmod')

    Vb = bond_potential_energy.bond_potential(molecule)
    Va = angle_potential_energy.angle_potential(molecule)
    Vd = dihedral_potential_energy.dihedral_potential(molecule)

    heuristic_conformation.heuristic_conformation(molecule, theta=np.pi/6, arquivo=path2+'_heuristic.mol2')

    Vb1 = bond_potential_energy.bond_potential(molecule)
    Va1 = angle_potential_energy.angle_potential(molecule)
    Vd1 = dihedral_potential_energy.dihedral_potential(molecule)
