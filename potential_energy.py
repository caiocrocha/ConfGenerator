import numpy as np
import matplotlib.pyplot as plt
from class_molecule import Molecule
from rotate import init_rotation_axis, rotate_atom, recursive_rotate_list

def calculate_pot_energy(molecule):
    Vb = 0
    natoms = molecule.num_atoms
    for atom in range(1, natoms):
        for i in range(atom):
            if(molecule.topology.bond_matrix[int(atom*(atom-1)/2 + i)]):
                btype = ('{}-{}'.format(molecule.topology.type[atom], 
                        molecule.topology.type[i])) # bond type (e.g. 'c3-c3' or 'c3-hc')
                try:
                    K = molecule.topology.bond_types[btype][0]  # elastic constant of the bond
                    b0 = molecule.topology.bond_types[btype][1] # equilibrium bond length
                except:
                    btype = ('{}-{}'.format(molecule.topology.type[i], 
                            molecule.topology.type[atom]))
                    try:
                        K = molecule.topology.bond_types[btype][0]
                        b0 = molecule.topology.bond_types[btype][1]
                    except: pass
                b = (((molecule.x[atom] - molecule.x[i])**2) + 
                     ((molecule.y[atom] - molecule.y[i])**2) + 
                     ((molecule.z[atom] - molecule.z[i])**2))**0.5
                d = b-b0
                Vb += K*(d**2) # acumulates the potential energy of the bond
    return Vb

def pot_energy_rotate(molecule, a, b, theta, ntimes, filename, pdf_name):
    molecule.write_pdb(filename, 'w+', 0)
    times = ntimes + 1
    degrees = np.zeros(times, dtype='float16')
    Vb = np.zeros(times, dtype='float32')
    degrees[0] = 0
    Vb[0] = calculate_pot_energy(molecule)
    rlist = []
    delta_theta = init_rotation_axis(molecule, a, b, theta, ntimes)
    recursive_rotate_list(molecule, b-1, a-1, rlist)
    for i in range(1, times):
        for atom in rlist:
            rotate_atom(molecule, atom)
        degrees[i] = delta_theta*i
        Vb[i] = calculate_pot_energy(molecule)
        molecule.write_pdb(filename, 'a', i)
    plt.plot(degrees, Vb)
    plt.xlabel('Degrees (rad)')
    plt.ylabel(r'Elastic potential (kcal $mol^{-1} \AA^{-2}$)')
    plt.grid()
    plt.suptitle('Rotation degrees X Elastic potential (Vb) graphic')
    plt.title('Number of rotations: {}'.format(ntimes), fontsize=10, loc='right')
    plt.savefig(pdf_name)
    plt.show()
    return Vb

molecule = Molecule()

#pdb, psf, out_pdb, out_psf = molecule.get_cmd_line()

molecule.read_psf("butane.psf")
molecule.read_pdb("butane_opt.pdb")
molecule.read_frcmod("butane.frcmod")
#Vb = calculate_pot_energy(molecule)
Vb = pot_energy_rotate(molecule, a=2, b=3, theta=2*np.pi, ntimes=360, 
                      filename='butane_opt1.pdb', pdf_name='Ep_butane_opt_360.pdf')
