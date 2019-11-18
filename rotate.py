import numpy as np
import argparse

class Molecule:
    def __init__(self):
        self.topology = self.Topology()
        self.id = list()    # "ATOM " or "HETATM"
        self.atom = list()  # atom name
        self.alt_loc = list()   # alternate location indicator
        self.residue = list()   # residue name
        self.chain = list() # chain identifier
        self.res_num = list()   # residue sequence number
        self.ins_code = list()   # insertion code for residues
        self.x = list() # orthogonal coordinates for X (in Angstroms)
        self.y = list() # orthogonal coordinates for Y (in Angstroms)
        self.z = list() # orthogonal coordinates for Z (in Angstroms)
        self.occup = list() # occupancy
        self.temp = list()  # temperature factor
        self.element = list()   # element symbol
        self.charge = list()

    def read_pdb(self, filename):
        with open(filename, "r") as inf:
            for line in inf:
                if(('ATOM' in line) or ('HETATM' in line)):
                    self.id.append(line[0:6].strip())
                    self.atom.append(line[12:16].strip())
                    self.alt_loc.append(line[16:17].strip())
                    self.residue.append(line[17:20].strip())
                    self.chain.append(line[21:22].strip())
                    self.res_num.append(int(line[22:26]))
                    self.ins_code.append(line[26:27].strip())
                    self.x.append(float(line[30:38]))
                    self.y.append(float(line[38:46]))
                    self.z.append(float(line[46:54]))
                    self.occup.append(float(line[54:60]))
                    self.temp.append(float(line[60:66].strip()))
                    self.element.append(line[72:75].strip())
                    self.charge.append(line[75:77].strip())
                elif(('TER' in line) or ('END' in line)):
                    break
                    
    def write_pdb(self, filename):
        with open(filename, "w+") as outf:
            outf.write("CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n")
            for i in range(len(self.id)):
                if(self.id[i] == 'ATOM' or self.id[i] == 'HETATM'):
                    outf.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}      {:>2s}{:2s}\n".format(
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
            outf.write("END\n")
    
    class Topology:
        def __init__(self):
            self.num_atoms = 0  # number of atoms
            self.segid = list()    # segment name
            self.resid = list()   # residue id
            self.resname = list()   # residue name
            self.name = list() # atom name
            self.type = list() # atom type
            self.charge = list() # charge
            self.mass = list()  # mass
            self.num_bonds = 0 # number of chemical bonds
            self.bond_list = list() # list of bonds - in pairs
            self.num_angles = 0   # number of angles
            self.angle_list = list()    # list of angles - in triplets
            self.num_dihedrals = 0  # number of dihedrals (torsions)
            self.dihedral_list = list() # list of dihedrals - in quartets
            self.bond_types = list()
            self.k_elastic = list() # constants of proportionality for the elastic potential energy equations 
            self.bond_lengths = list() # bond lengths
            self.angles = list()    # angle values
    
        def read_psf(self, filename):
            with open(filename, "r") as inf:
                line = inf.readline()
                while line:
                    if("!NATOM" in line):
                        self.num_atoms = int(line.split()[0])
                        for i in range(self.num_atoms):
                            line = inf.readline().split()
                            self.segid.append(line[1])
                            self.resid.append(int(line[2]))
                            self.resname.append(line[3])
                            self.name.append(line[4])
                            self.type.append(line[5])
                            self.charge.append(float(line[6]))
                            self.mass.append(float(line[7]))
                    elif("!NBOND" in line):
                        self.num_bonds = int(line.split()[0])
                        num_lines = int(self.num_bonds/4)
                        if(self.num_bonds%2 != 0):
                            num_lines+=1
                        for i in range(num_lines):
                            self.bond_list.extend(inf.readline().split())
                    elif("!NTHETA" in line):
                        self.num_angles = int(line.split()[0])
                        num_lines = int(self.num_angles/3)
                        if(self.num_angles%2 != 0):
                            num_lines+=1
                        for i in range(num_lines):
                            self.angle_list.extend(inf.readline().split())
                    elif("!NPHI" in line):
                        self.num_dihedrals = int(line.split()[0])
                        num_lines = int(self.num_dihedrals/2)
                        if(self.num_dihedrals%2 != 0):
                            num_lines+=1
                        for i in range(num_lines):
                            self.dihedral_list.extend(inf.readline().split())
                    line = inf.readline()
    
        def write_psf(self, filename):
            with open(filename, "w+") as outf:
                outf.write("PSF CHEQ EXT XPLOR\n\n{:10d} !NTITLE\n\n\n{:10d} !NATOM\n".format(1, self.num_atoms))
                for i in range(self.num_atoms):
                    outf.write("{:10d} {:5s}{:5d} {:>10s}     {:^5s}     {:5s}{:12.6f}  {:12.4f}\n".format(
                    i+1,
                    self.segid[i],
                    self.resid[i],
                    self.resname[i],
                    self.name[i],
                    self.type[i],
                    self.charge[i],
                    self.mass[i]))
                outf.write("\n{:10d} !NBOND: bonds\n".format(self.num_bonds))
                for i in range(self.num_bonds*2):
                    outf.write("{:>10s}".format(self.bond_list[i]))
                    if(((i+1)%8 == 0 and i != 0) or (i == self.num_bonds*2-1)):
                        outf.write("\n")
                outf.write("\n")
                outf.write("{:10d} !NTHETA: angles\n".format(self.num_angles))
                for i in range(self.num_angles*3):
                    outf.write("{:>10s}".format(self.angle_list[i]))
                    if(((i+1)%9 == 0 and i != 0) or (i == self.num_angles*3-1)):
                        outf.write("\n")
                outf.write("\n")
                outf.write("{:10d} !NPHI: dihedrals\n".format(self.num_dihedrals))
                for i in range(self.num_dihedrals*4):
                    outf.write("{:>10s}".format(self.dihedral_list[i]))
                    if(((i+1)%8 == 0 and i != 0) or (i == self.num_dihedrals*4-1)):
                        outf.write("\n")
                outf.write("\n")
                
        def read_frcmod(self, filename):
            with open(filename, "r") as inf:
                line = inf.readline()
                while line:
                    if "BOND" in line:
                        line = inf.readline()
                        while line.strip():
                            self.bond_types.append(line[0:5].strip())
                            self.k_elastic.append(float(line[5:13].strip()))
                            self.bond_lengths.append(float(line[13:21].strip()))
                            line = inf.readline()
                    line = inf.readline()
                
    def get_bond_list(self, a):
        blist = list()
        for i in range(0, self.topology.num_bonds*2, 2):
            if(str(a) in self.topology.bond_list[i:i+2]):
                blist.append(int(self.topology.bond_list[i]))
                blist.append(int(self.topology.bond_list[i+1]))
        return blist
    
    def get_angle_list(self, a):
        alist = list()
        for i in range(0, self.topology.num_angles*3, 3):
            if(a == int(self.topology.angle_list[i+1])):
                alist.append(int(self.topology.angle_list[i]))
                alist.append(int(self.topology.angle_list[i+1]))
                alist.append(int(self.topology.angle_list[i+2]))
        return alist
    
    def get_dihedral_list(self, a):
        dlist = list()
        for i in range(0, self.topology.num_dihedrals*4, 4):
            if(str(a) in self.topology.dihedral_list[i+1:i+3]):
                dlist.append(int(self.topology.dihedral_list[i]))
                dlist.append(int(self.topology.dihedral_list[i+1]))
                dlist.append(int(self.topology.dihedral_list[i+2]))
                dlist.append(int(self.topology.dihedral_list[i+3]))
        return dlist
    
    def rotate_list(self, r, c, start, b):
        if(c < len(r) < self.topology.num_atoms):
            for i in range(start, self.topology.num_bonds*2, 2):
                p = int(self.topology.bond_list[i])
                q = int(self.topology.bond_list[i+1])
                e = r[c]
                if(e == p and q != b):
                    r.append(q)
                    start = i+2
                    break
                elif(e == q and p != b): 
                    r.append(p)
                    start = i+2
                    break
            for i in range(start, self.topology.num_bonds*2, 2):
                p = int(self.topology.bond_list[i])
                q = int(self.topology.bond_list[i+1])
                e = r[c]
                if(e == p and q != b):
                    r.append(q)
                elif(e == q and p != b): 
                    r.append(p)
            self.rotate_list(r, c+1, start, b)
    
    def rotate(self, a, b, theta):   # bond a-b, rotate chain from b // theta in radians
        r = list()  # list of atoms to be rotated
        for i in range(0, self.topology.num_bonds*2, 2):
            p = int(self.topology.bond_list[i])
            q = int(self.topology.bond_list[i+1])
            if(p == b and q != a):
                r.append(q)
                start = i+2
                break
            elif(q == b and p != a): 
                r.append(p)
                start = i+2
                break
        for i in range(start, self.topology.num_bonds*2, 2):
            p = int(self.topology.bond_list[i])
            q = int(self.topology.bond_list[i+1])
            if(p == b and q != a):
                r.append(q)
            elif(q == b and p != a): 
                r.append(p)
        self.rotate_list(r, 0, start, b)
        r = list(set(r)) # removes duplicates
        '''
        For the construction of the rotation matrix around the axis defined by the bond that
        connects atoms a and b (vector ab, a to b), there needs to be determined ab's unit vector, 
        which is mathematicaly defined by its original components divided by the square root
        of its magnitude.
        '''
        Ax = self.x[a-1]
        Ay = self.y[a-1]
        Az = self.z[a-1]
        Bx = self.x[b-1]
        By = self.y[b-1]
        Bz = self.z[b-1]
        M = (((Bx - Ax)**2) + ((By - Ay)**2) + ((Bz - Az)**2))**0.5
        # M is the magnitude of the unit vector u
        ux = (Bx - Ax)/M # x component of u
        uy = (By - Ay)/M # y component of u
        uz = (Bz - Az)/M # z component of u
        cos0 = np.cos(theta) # cosine of the angle of rotation
        sin0 = np.sin(theta) # sine of the angle of rotation
        t = 1 - cos0
        for i in range(len(r)):
            j = r[i]-1
            x1 = self.x[j] - Ax
            y1 = self.y[j] - Ay
            z1 = self.z[j] - Az
            self.x[j] = round(((cos0 + t*(ux**2))*x1 + (t*ux*uy - sin0*uz)*y1 + (t*ux*uz + sin0*uy)*z1), 3) + Ax
            self.y[j] = round(((t*ux*uy + sin0*uz)*x1 + (cos0 + t*(uy**2))*y1 + (t*uy*uz - sin0*ux)*z1), 3) + Ay
            self.z[j] = round(((t*ux*uz - sin0*uy)*x1 + (t*uy*uz + sin0*ux)*y1 + (cos0 + t*(uz**2))*z1), 3) + Az
            # rotation matrix applied to point
            
    def calculate_pot_energy(self):
        Vb = 0  # total elastic potential energy in the molecule
        for i in range(0, self.topology.num_bonds*2, 2):
            for j in range(len(self.topology.bond_types)):
                bond = (self.topology.bond_types[j].split('-')) # bond type (e.g. 'c3-c3' or 'c3-hc')
                ind_a = int(self.topology.bond_list[i])-1 # index of atom a
                ind_b = int(self.topology.bond_list[i+1])-1 # index of atom b
                atom_a = self.topology.type[ind_a] # atom type (e.g. 'c3' or 'hc')
                atom_b = self.topology.type[ind_b]
                if((atom_a == bond[0] and atom_b == bond[1]) or (atom_a == bond[1] and atom_b == bond[0])):
                    K = self.topology.k_elastic[j] # elastic constant of the bond
                    b = (((self.x[ind_a]-self.x[ind_b])**2)+((self.y[ind_a]-self.y[ind_b])**2)+((self.z[ind_a]-self.z[ind_b])**2))**0.5
                    # b is the distance between atoms a and b
                    b0 = self.topology.bond_lengths[j] # b0 is the equilibrium bond length
                    Vb += K*((b-b0)**2) # acumulates the potential energy of the bond
        return Vb

'''
def get_cmd_line():
    parser = argparse.ArgumentParser(description = 'PDB and PSF reader, chain rotator.')

    parser.add_argument('--pdb', action = 'store',     dest = 'pdb', required = True, help = 'PDB file.')
    parser.add_argument('--psf', action = 'store',     dest = 'psf', required = True, help = 'PSF file.')
    parser.add_argument('--new_pdb', action = 'store',     dest = 'new_pdb', required = True, help = 'New PDB file.')
    parser.add_argument('--new_psf', action = 'store',     dest = 'new_psf', required = True, help = 'New PSF file.')

    arguments = parser.parse_args()


    # Assign from cmd line.
    return (arguments.pdb,
            arguments.psf,
            arguments.new_pdb,
            arguments.new_psf)
'''           
        
molecule = Molecule()

# pdb, psf, out_pdb, out_psf = get_cmd_line()

molecule.read_pdb("butane_opt.pdb")
molecule.topology.read_psf("butane.psf")
molecule.topology.read_frcmod("butane.frcmod")
bond_list_C3 = molecule.get_bond_list(3)
angle_list_C3 = molecule.get_angle_list(3)
dihed_list_C3 = molecule.get_dihedral_list(3)
molecule.rotate(2, 3, np.pi/2)
Vb = molecule.calculate_pot_energy()
print("Vb =", Vb)
#molecule.write_pdb("butane_opt_rotated.pdb")

'''
for i in range(0, len(bond_list_C3), 2):
	print(bond_list_C3[i:i+2])
print("\n")
for i in range(0, len(angle_list_C3), 3):
	print(angle_list_C3[i:i+3])
print("\n")
for i in range(0, len(dihed_list_C3), 4):
	print(dihed_list_C3[i:i+4])
print("\n")
for i in range(0, len(molecule.topology.bonds), 2):
	print(molecule.topology.bond_list[i:i+2])
print("\n")
'''