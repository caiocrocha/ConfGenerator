import argparse
import numpy as np
import math
from class_trimatrix import Trimatrix

def is_main():
    return __name__ == '__main__'

class Molecule:
    def __init__(self):
        self.topology = self.Topology()
        self.id        = [] # 'ATOM' or 'HETATM'
        self.atom      = [] # atom name
        self.alt_loc   = [] # alternate location indicator
        self.residue   = [] # residue name
        self.chain     = [] # chain identifier
        self.res_num   = [] # residue sequence number
        self.ins_code  = [] # insertion code for residues
        self.x         = [] # orthogonal coordinates for X (in Angstroms)
        self.y         = [] # orthogonal coordinates for Y (in Angstroms)
        self.z         = [] # orthogonal coordinates for Z (in Angstroms)
        self.occup     = [] # occupancy
        self.temp      = [] # temperature factor
        self.element   = [] # element symbol
        self.charge    = []
        self.num_atoms = 0  # number of atoms
        
    class Topology:
        def __init__(self):
            self.segid       = [] # segment name
            self.resid       = [] # residue id
            self.resname     = [] # residue name
            self.name        = [] # atom name
            self.type        = [] # atom type
            self.charge      = [] # charge
            self.mass        = [] # mass
            self.num_bonds   = 0  # number of chemical bonds
            self.bond_list   = []
            self.bond_matrix = None
            # Matrix of bonds, which is a symmetric matrix in its linear 
            # representation, for the sake of memory savings. The presence or 
            # absence of a bond between atoms (row, colune) is indicated in 
            # self.bond_matrix[k], where k corresponds to indices (i, j) in the 
            # matricial representation.

            self.num_angles     = 0  # number of angles
            self.angle_list     = [] # list of angles - in triplets
            self.num_dihedrals  = 0  # number of dihedrals (torsions)
            self.dihedral_list  = [] # list of dihedrals - in quartets
            self.bond_types     = {}
            # e.g. self.bond_types['c3-hc'][0] or self.bond_types['c3-hc'][1], where 
            # index 0 refers to the "spring" constant for the elastic potential energy 
            # equations and index 1 refers to the equilibrium bond length.
            self.angle_types    = {}
            # e.g. self.angle_types['hc-c3-hc'][0] or self.bond_types['hc-c3-hc'][1], 
            # where index 0 refers to the angle "spring" constant and index 1 refers 
            # to the equilibrium angle converted to radians.
            self.dihedral_types = {}
            # e.g. self.dihedral_types['c3-c3-c3-c3'][0], 
            # self.dihedral_types['c3-c3-c3-c3'][1] (dihedral angle "spring" constant), 
            # self.dihedral_types['c3-c3-c3-c3'][2] (dihedral angle converted to radians) 
            # or self.dihedral_types['c3-c3-c3-c3'][3] (multiplicity).

    # def get_cmd_line():
    #     parser = argparse.ArgumentParser(description = 'PDB and PSF reader, chain rotator.')
    
    #     parser.add_argument('--pdb', action = 'store',     dest = 'pdb', required = True, help = 'PDB file.')
    #     parser.add_argument('--psf', action = 'store',     dest = 'psf', required = True, help = 'PSF file.')
    #     parser.add_argument('--new_pdb', action = 'store',     dest = 'new_pdb', required = True, help = 'New PDB file.')
    #     parser.add_argument('--new_psf', action = 'store',     dest = 'new_psf', required = True, help = 'New PSF file.')
    
    #     arguments = parser.parse_args()
    
    #     # Assign from cmd line.
    #     return (arguments.pdb,
    #             arguments.psf,
    #             arguments.new_pdb,
    #             arguments.new_psf)  

    def read_mol2(self, filename):
        with open (filename, 'r') as inf:
            line = inf.readline()
            while line:
                if '@<TRIPOS>MOLECULE' in line:
                    line = inf.readline()
                    line = inf.readline().split()
                    self.num_atoms = int(line[0])
                    self.topology.num_bonds = int(line[1])
                elif '@<TRIPOS>ATOM' in line:
                    natoms = self.num_atoms
                    self.atom = np.empty(natoms, dtype='U4')
                    self.x = np.zeros(natoms, dtype='float32')
                    self.y = np.zeros(natoms, dtype='float32')
                    self.z = np.zeros(natoms, dtype='float32')
                    self.topology.type = np.empty(natoms, dtype='U5')
                    self.topology.charge = np.zeros(natoms, dtype='float16')
                    for i in range(natoms):
                        line = inf.readline().split()
                        self.atom[i] = line[1]
                        self.x[i] = float(line[2])
                        self.y[i] = float(line[3])
                        self.z[i] = float(line[4])
                        self.topology.type[i] = line[5]
                        self.topology.charge[i] = line[8]
                elif '@<TRIPOS>BOND' in line:
                    nbonds = self.topology.num_bonds
                    self.topology.bond_list = np.zeros(nbonds * 2, dtype='uint8')
                    self.topology.bond_matrix = np.zeros(Trimatrix.get_size(
                        self.num_atoms
                    ), dtype=bool)
                    for i in range(0, nbonds*2, 2):
                        line = inf.readline().split()
                        self.topology.bond_matrix[Trimatrix.get_index(
                            int(line[1]) - 1, int(line[2]) - 1
                        )] = True
                        self.topology.bond_list[i] = line[1]
                        self.topology.bond_list[i+1] = line[2]
                line = inf.readline()

    def write_mol2(self, filename):
        with open(filename, 'w+') as outf:
            natoms = self.num_atoms
            nbonds = self.topology.num_bonds
            outf.write('@<TRIPOS>MOLECULE\nMOL\n')
            outf.write('{:>5d}{:>5d}{:>5}{:>5}{:>5}\n'.format(natoms, nbonds, '', '', ''))
            outf.write('{}\n{}\n\n\n'.format('', ''))
            outf.write('@<TRIPOS>ATOM\n')
            for i in range(natoms):
                outf.write('{0:>4} {1:>4} {2:>13.4f} {3:>9.4f} {4:>9.4f} {5:>4} {6} {7} {8:>7.4f}\n'.format(
                    i+1, self.atom[i], self.x[i], self.y[i], self.z[i], self.topology.type[i], '', '', self.topology.charge[i]
                ))
            outf.write('@<TRIPOS>BOND\n')
            count = 0
            for i in range(nbonds):
                outf.write('{0:>5} {1:>5} {2:>5} {3:>2}\n'.format(
                    i+1, self.topology.bond_list[count], self.topology.bond_list[count+1], 1
                ))
                count += 2

    def read_psf(self, filename):
        with open(filename, 'r') as inf:
            line = inf.readline()
            while line:
                if'!NATOM' in line:
                    natoms = self.num_atoms = int(line.split()[0])
                    self.topology.segid     = np.empty(natoms, dtype='U3')
                    self.topology.resid     = np.zeros(natoms, dtype='uint8')
                    self.topology.resname   = np.empty(natoms, dtype='U3')
                    self.topology.name      = np.empty(natoms, dtype='U5')
                    self.topology.type      = np.empty(natoms, dtype='U5')
                    self.topology.charge    = np.zeros(natoms, dtype='float16')
                    self.topology.mass      = np.zeros(natoms, dtype='float16')

                    for i in range(self.num_atoms):
                        line = inf.readline().split()
                        self.topology.segid[i]   = line[1]
                        self.topology.resid[i]   = line[2]
                        self.topology.resname[i] = line[3]
                        self.topology.name[i]    = line[4]
                        self.topology.type[i]    = line[5]
                        self.topology.charge[i]  = line[6]
                        self.topology.mass[i]    = line[7]

                elif('!NBOND' in line):
                    nbonds = int(line.split()[0])
                    self.topology.num_bonds = nbonds
                    num_lines = math.ceil(nbonds/4)
                    self.topology.bond_list   = np.zeros(nbonds*2, dtype='uint8')
                    self.topology.bond_matrix = np.zeros(Trimatrix.get_size(
                        self.num_atoms
                        ), dtype=bool)
                    count = 0

                    for i in range(num_lines):
                        line = inf.readline().split()
                        for j in range(0, len(line), 2):
                            self.topology.bond_matrix[Trimatrix.get_index(
                                int(line[j])-1, int(line[j+1])-1
                                )] = True
                            self.topology.bond_list[count] = line[j]
                            self.topology.bond_list[count+1] = line[j+1]
                            count += 2

                elif('!NTHETA' in line):
                    nangles = self.topology.num_angles = int(line.split()[0])
                    num_lines = math.ceil(nangles/3)
                    self.topology.angle_list = np.zeros(nangles*3, dtype='uint8')
                    count = 0
                    for i in range(num_lines):
                        line = inf.readline().split()
                        for j in range(len(line)):
                            self.topology.angle_list[count] = line[j]
                            count+=1

                elif('!NPHI' in line):
                    ndihed = self.topology.num_dihedrals = int(line.split()[0])
                    num_lines = math.ceil(ndihed/2)
                    self.topology.dihedral_list = np.zeros(ndihed*4, dtype='uint8')
                    count = 0
                    for i in range(num_lines):
                        line = inf.readline().split()
                        for j in range(len(line)):
                            self.topology.dihedral_list[count] = line[j]
                            count+=1
                line = inf.readline()

    def write_psf(self, filename):
        with open(filename, 'w+') as outf:
            outf.write('PSF CHEQ EXT XPLOR\n\n{:10d} !NTITLE\n\n\n'.format(1))
            natoms = self.num_atoms
            outf.write('{:10d} !NATOM\n'.format(natoms))

            for i in range(self.num_atoms):
                outf.write('{:10d} {:5s}{:5d} {:>10s}     {:^5s}     {:5s}{:12.6f}  {:12.4f}\n'.format(
                i+1,
                self.topology.segid[i],
                self.topology.resid[i],
                self.topology.resname[i],
                self.topology.name[i],
                self.topology.type[i],
                self.topology.charge[i],
                self.topology.mass[i])
                    )
            nbonds = self.topology.num_bonds
            outf.write('\n{:10d} !NBOND: bonds\n'.format(nbonds))

            for i in range(nbonds*2):
                outf.write('{:>10d}'.format(self.topology.bond_list[i]))
                if ((i + 1) % 8 == 0 and i != 0) or (i == nbonds * 2 - 1):
                    outf.write('\n')

            outf.write('\n')
            nangles = self.topology.num_angles
            outf.write('{:10d} !NTHETA: angles\n'.format(nangles))

            for i in range(nangles*3):
                outf.write('{:>10d}'.format(self.topology.angle_list[i]))
                if ((i + 1) % 9 == 0 and i != 0) or (i == nangles * 3 - 1):
                    outf.write('\n')

            outf.write('\n')
            ndihed = self.topology.num_dihedrals
            outf.write('{:10d} !NPHI: dihedrals\n'.format(ndihed))

            for i in range(ndihed*4):
                outf.write('{:>10d}'.format(self.topology.dihedral_list[i]))
                if ((i + 1) % 8 == 0 and i != 0) or (i == ndihed * 4 - 1):
                    outf.write('\n')

            outf.write('\n')
            outf.write('{:10d} !NIMPHI: impropers\n\n\n'.format(0))
            outf.write('{:10d} !NDON: donors\n\n\n'.format(0))
            outf.write('{:10d} !NACC: acceptors\n\n\n'.format(0))
            outf.write('{:10d} !NNB\n\n'.format(0))

            for i in range(natoms):
               outf.write('{:>10d}'.format(0))
               if ((i + 1) % 8 == 0 and i != 0) or (i == natoms - 1):
                    outf.write('\n')
            
    def read_pdb(self, filename):
        with open(filename, 'r') as inf:
            if self.num_atoms != 0:
                # Statically allocates memory as long as there has already been a 
                # corresponding .psf file read, which will determine to the function 
                # the number of spaces (number of atoms) that need to be preallocated. 
                # Static allocation is preferable over dynamic allocation, since it
                # uses numpy arrays, which can save up memory by limiting the variables 
                # sizes. Also, the 'append' method is slower than preallocating the 
                # elements and then assigning them values.
                natoms        = self.num_atoms
                self.id       = np.empty(natoms, dtype='U6')
                self.atom     = np.empty(natoms, dtype='U4')
                self.alt_loc  = np.empty(natoms, dtype='U1')
                self.residue  = np.empty(natoms, dtype='U3')
                self.chain    = np.empty(natoms, dtype='U1')
                self.res_num  = np.zeros(natoms, dtype='uint8')
                self.ins_code = np.empty(natoms, dtype='U1')
                self.x        = np.zeros(natoms, dtype='float32')
                self.y        = np.zeros(natoms, dtype='float32')
                self.z        = np.zeros(natoms, dtype='float32')
                self.occup    = np.zeros(natoms, dtype='float16')
                self.temp     = np.zeros(natoms, dtype='float16')
                self.element  = np.empty(natoms, dtype='U3')
                self.charge   = np.empty(natoms, dtype='U2')
                i = 0

                for line in inf:
                    if 'ATOM' in line or 'HETATM' in line:
                        self.id[i]       = line[0:6].strip()
                        self.atom[i]     = line[12:16].strip()
                        self.alt_loc[i]  = line[16:17].strip()
                        self.residue[i]  = line[17:20].strip()
                        self.chain[i]    = line[21:22].strip()
                        self.res_num[i]  = line[22:26]
                        self.ins_code[i] = line[26:27].strip()
                        self.x[i]        = line[30:38]
                        self.y[i]        = line[38:46]
                        self.z[i]        = line[46:54]
                        self.occup[i]    = line[54:60]
                        self.temp[i]     = line[60:66].strip()
                        self.element[i]  = line[72:75].strip()
                        self.charge[i]   = line[75:77].strip()
                        i+=1

                    elif 'TER' in line or 'END' in line:
                        break

            else:
                for line in inf:
                    if 'ATOM' in line or 'HETATM' in line:
                        self.id.append(         line[0:6].strip())
                        self.atom.append(       line[12:16].strip())
                        self.alt_loc.append(    line[16:17].strip())
                        self.residue.append(    line[17:20].strip())
                        self.chain.append(      line[21:22].strip())
                        self.res_num.append(int(line[22:26]))
                        self.ins_code.append(   line[26:27].strip())
                        self.x.append(float(    line[30:38]))
                        self.y.append(float(    line[38:46]))
                        self.z.append(float(    line[46:54]))
                        self.occup.append(float(line[54:60]))
                        self.temp.append(float( line[60:66].strip()))
                        self.element.append(    line[72:75].strip())
                        self.charge.append(     line[75:77].strip())
                    elif 'TER' in line or 'END' in line:
                        break
                self.num_atoms = len(self.id)
                    
    def write_pdb(self, filename, mode, n):
        with open(filename, mode) as outf:
            outf.write('MODEL     {:4d}\n'.format(n))
            # outf.write('CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n')
            for i in range(self.num_atoms):
                if(self.id[i] == 'ATOM' or self.id[i] == 'HETATM'):
                    outf.write('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:>2s}{:2s}\n'.format(
                    self.id[i],
                    i+1,    # atom serial number
                    self.atom[i],
                    self.alt_loc[i],
                    self.residue[i],
                    self.chain[i],
                    self.res_num[i],
                    self.ins_code[i],
                    self.x[i],
                    self.y[i],
                    self.z[i],
                    self.occup[i],
                    self.temp[i],
                    self.element[i],
                    self.charge[i]))
            outf.write('END\n')

    def read_frcmod(self, filename):
        with open(filename, 'r') as inf:
            line = inf.readline()
            while line:
                if 'BOND' in line:
                    line = inf.readline()
                    while line.strip():
                        self.topology.bond_types[line[0:5].strip()] = (
                            (float(line[5:13].strip()), float(line[13:21].strip()))
                            )
                        line = inf.readline()
                elif 'ANGLE' in line:
                    line = inf.readline()
                    while line.strip():
                        angle_rad = float(line[17:29].strip())*np.pi/180
                        self.topology.angle_types[line[0:8].strip()] = (
                            (float(line[8:17].strip()), angle_rad)
                            )
                        line = inf.readline()
                elif 'DIHE' in line:
                    line = inf.readline()
                    while line.strip():
                        key = line[0:11]
                        angle_rad = float(line[24:38].strip())*np.pi/180
                        line = [int(line[11:15].strip()), float(line[15:24].strip()), 
                                 angle_rad, float(line[38:54].strip())]
                        if key not in self.topology.dihedral_types.keys():
                            self.topology.dihedral_types[key] = [line]
                        else:
                            self.topology.dihedral_types[key].append(line)
                        line = inf.readline()
                line = inf.readline()
                
    def get_bond_list(self, a):
        blist = []
        a = a-1
        k = int(a*(a-1)/2)

        for i in range(a):
            if(self.topology.bond_matrix[k+i]):
                blist.append(i+1)

        for i in range(a+1, self.num_atoms):
            if(self.topology.bond_matrix[int(i*(i-1)/2 + a)]):
                blist.append(i+1)

        return blist
    
    def get_angle_list(self, a):
        alist = []

        for i in range(0, self.topology.num_angles*3, 3):
            if(a == int(self.topology.angle_list[i+1])):
                alist.append(int(self.topology.angle_list[i]))
                alist.append(int(self.topology.angle_list[i+1]))
                alist.append(int(self.topology.angle_list[i+2]))

        return alist
    
    def get_dihedral_list(self, a):
        dlist = []

        for i in range(0, self.topology.num_dihedrals*4, 4):
            if(a in self.topology.dihedral_list[i+1:i+3]):
                dlist.append(int(self.topology.dihedral_list[i]))
                dlist.append(int(self.topology.dihedral_list[i+1]))
                dlist.append(int(self.topology.dihedral_list[i+2]))
                dlist.append(int(self.topology.dihedral_list[i+3]))

        return dlist
    
##############################################################################
if is_main():
    # code here
    pass
