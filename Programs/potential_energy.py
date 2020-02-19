import argparse
import numpy as np
import os

from Classes.Molecule import Molecule
import DihedralPotential.dihedral_potential_energy as dihedral_potential_energy
import AnglePotential.angle_potential_energy as angle_potential_energy
import BondPotential.bond_potential_energy as bond_potential_energy
import BondPotential.minimize_bond_energy as minimize_bond_energy
import BondForce.bond_force as bond_force
import AngleForce.angle_force as angle_force
import DihedralForce.dihedral_force as dihedral_force
import DihedralPotential.dihedral_energy_graphic as dihedral_energy_graphic
import heuristic_conformation

def is_main():
    return __name__ == '__main__'

##############################################################################
if is_main():

    molecule = Molecule()

	path = os.getcwd()

	'''
    path1 = path + '../butane'
    path2 = path1 + '/sqm/sqm'
    path3 = path1 + '/butane'
    '''
    
    path1 = path + '../3-metil-pentano'
    path2 = path1 + '/3-metil-pentano'
    path3 = path1 + '/3-metil-pentano'

    molecule.read_mol2(path2 + '.mol2')
    molecule.gen_dihed_list_from_angle_list()
    molecule.read_frcmod(path3 + '.frcmod')
    
    Vb = bond_potential_energy.bond_potential(molecule)
    Va = angle_potential_energy.angle_potential(molecule)
    Vd = dihedral_potential_energy.total_dihedral_potential(molecule)

    heuristic_conformation.heuristic_conformation(molecule, theta=np.pi/6)
    molecule.write_mol2(path2+'_heuristic.mol2')

    Vb1 = bond_potential_energy.bond_potential(molecule)
    Va1 = angle_potential_energy.angle_potential(molecule)
    Vd1 = dihedral_potential_energy.total_dihedral_potential(molecule)
