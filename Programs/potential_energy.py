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

def get_cmd_line():
    parser = argparse.ArgumentParser(description='MOL2 and FRCMOD reader, potential energy calculator.')
    parser.add_argument('--mol2', 	  action='store', dest='mol2',	   required=True, help='MOL2 file.')
    parser.add_argument('--frcmod',   action='store', dest='frcmod',   required=True, help='FRCMOD file.')
    parser.add_argument('--new_mol2', action='store', dest='new_mol2', required=True, help='New MOL2 file.')
    arguments = parser.parse_args()
    # Assign from cmd line.
    return arguments.mol2, arguments.frcmod, arguments.new_mol2

def main():

    molecule = Molecule()

    path = os.path.dirname(os.path.abspath(__file__))

    # mol2, frcmod, new_mol2 = get_cmd_line()

    '''
    path1 = path + '/../3-metil-pentano'
    path2 = path3 = path1 + '/3-metil-pentano'

    path1 = path + '/../butane'
    path2 = path1 + '/sqm/sqm'
    path3 = path1 + '/butane_parm'
    '''

    path1 = path + '/../Etileno_glicol'
    path2 = path3 = path1 + '/ligand'

    molecule.read_mol2(path2 + '.mol2')
    molecule.gen_dihed_list_from_angle_list()
    molecule.read_frcmod(path3 + '.frcmod')
    '''
    Vb = bond_potential_energy.bond_potential(molecule)
    Va = angle_potential_energy.angle_potential(molecule)
    Vd = dihedral_potential_energy.total_dihedral_potential(molecule)

    heuristic_conformation.heuristic_conformation(molecule, np.pi/6)

    Vd1 = dihedral_potential_energy.total_dihedral_potential(molecule)
    print('Vb =', Vb)
    print('Va =', Va)
    print('Vd (before) =', Vd)
    print('Vd (after) =', Vd1)
    '''
    dihedral_energy_graphic.dihedral_energy_graphic(molecule, 4, 1, 2, 3, np.pi/180, 360,
		pdf_name=path+'/../Graphs/PE_ligand_semH.pdf', write_mol2=False, mol2_name=path2 + '_rotated.mol2')

if __name__ == '__main__': main()
