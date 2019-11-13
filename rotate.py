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

    def readPDB(self, filename):
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
                    
    def writePDB(self, filename):
        with open(filename, "w+") as outf:
            #outf.write("CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1\n")
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
            self.bonds = list() # list of bonds - in pairs
            self.num_angles = 0   # number of angles
            self.angles = list()    # list of angles - in triplets
            self.num_dihedrals = 0  # number of dihedrals(torsions)
            self.dihedrals = list() # list of dihedrals - in quartets
    
        def readPSF(self, filename):
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
                            self.bonds.extend(inf.readline().split())
                    elif("!NTHETA" in line):
                        self.num_angles = int(line.split()[0])
                        num_lines = int(self.num_angles/3)
                        if(self.num_angles%2 != 0):
                            num_lines+=1
                        for i in range(num_lines):
                            self.angles.extend(inf.readline().split())
                    elif("!NPHI" in line):
                        self.num_dihedrals = int(line.split()[0])
                        num_lines = int(self.num_dihedrals/2)
                        if(self.num_dihedrals%2 != 0):
                            num_lines+=1
                        for i in range(num_lines):
                            self.dihedrals.extend(inf.readline().split())
                    line = inf.readline()
    
        def writePSF(self, filename):
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
                    outf.write("{:>10s}".format(self.bonds[i]))
                    if(((i+1)%8 == 0 and i != 0) or (i == self.num_bonds*2-1)):
                        outf.write("\n")
                outf.write("\n")
                outf.write("{:10d} !NTHETA: angles\n".format(self.num_angles))
                for i in range(self.num_angles*3):
                    outf.write("{:>10s}".format(self.angles[i]))
                    if(((i+1)%9 == 0 and i != 0) or (i == self.num_angles*3-1)):
                        outf.write("\n")
                outf.write("\n")
                outf.write("{:10d} !NPHI: dihedrals\n".format(self.num_dihedrals))
                for i in range(self.num_dihedrals*4):
                    outf.write("{:>10s}".format(self.dihedrals[i]))
                    if(((i+1)%8 == 0 and i != 0) or (i == self.num_dihedrals*4-1)):
                        outf.write("\n")
                outf.write("\n")
                
    def bond_list(self, a):
        bond_list = list()
        for i in range(0, len(self.topology.bonds), 2):
            if(a == int(self.topology.bonds[i]) or a == int(self.topology.bonds[i+1])):
                bond_list.append(int(self.topology.bonds[i]))
                bond_list.append(int(self.topology.bonds[i+1]))
        return bond_list
    
    def angle_list(self, a):
        angle_list = list()
        for i in range(0, len(self.topology.angles), 3):
            if(a == int(self.topology.angles[i+1])):
                angle_list.append(int(self.topology.angles[i]))
                angle_list.append(int(self.topology.angles[i+1]))
                angle_list.append(int(self.topology.angles[i+2]))
        return angle_list
    
    def dihedral_list(self, a):
        dihedral_list = list()
        for i in range(0, len(self.topology.dihedrals), 4):
            if(a == int(self.topology.dihedrals[i+1]) or a == int(self.topology.dihedrals[i+2])):
                dihedral_list.append(int(self.topology.dihedrals[i]))
                dihedral_list.append(int(self.topology.dihedrals[i+1]))
                dihedral_list.append(int(self.topology.dihedrals[i+2]))
                dihedral_list.append(int(self.topology.dihedrals[i+3]))
        return dihedral_list

    def rotate(self, a, b, theta):   # bond a-b, rotate chain from b // theta in radians
        r = list()  # list of atoms to be rotated
        for i in range(0, len(self.topology.bonds), 2):
            if(b == int(self.topology.bonds[i]) and a != int(self.topology.bonds[i+1])):
                r.append(int(self.topology.bonds[i+1]))
                start = i+2
                break
            elif(b == int(self.topology.bonds[i+1]) and a != int(self.topology.bonds[i])): 
                r.append(int(self.topology.bonds[i]))
                start = i+2
                break
        for i in range(start, len(self.topology.bonds), 2):
            if(b == int(self.topology.bonds[i]) and a != int(self.topology.bonds[i+1])):
                r.append(int(self.topology.bonds[i+1]))
            elif(b == int(self.topology.bonds[i+1]) and a != int(self.topology.bonds[i])): 
                r.append(int(self.topology.bonds[i]))
        self.rotate_list(r, 0, start, b)
        r = list(set(r))
        '''
        For the construction of the rotation matrix around the axis defined by the bond that
        connects atoms a and b (vector ab, a to b), there needs to be determined ab's unit vector, 
        which is mathematicaly defined by its original components divided by the square root
        of its magnitude.
        '''
        m = (((self.x[b-1] - self.x[a-1])**2) + ((self.y[b-1] - self.y[a-1])**2) + ((self.z[b-1] - self.z[a-1])**2))**0.5
        # m is the magnitude of the unit vector u
        ux = (self.x[b-1] - self.x[a-1])/m # x component of u
        uy = (self.y[b-1] - self.y[a-1])/m # y component of u
        uz = (self.z[b-1] - self.z[a-1])/m # z component of u
        cos0 = np.cos(theta) # cosine of the angle of rotation
        sin0 = np.sin(theta) # sine of the angle of rotation
        for i in range(len(r)):
            x1 = self.x[r[i]-1]
            y1 = self.y[r[i]-1]
            z1 = self.z[r[i]-1]
            self.x[r[i]-1] = round(((cos0 + (1-cos0)*(ux**2))*x1 + ((1-cos0)*ux*uy + sin0*uz)*y1 + ((1-cos0)*ux*uz + sin0*uy)*z1), 3)
            self.y[r[i]-1] = round((((1-cos0)*ux*uy + sin0*uz)*x1 + (cos0 + (1-cos0)*(uy**2))*y1 + ((1-cos0)*uy*uz + sin0*ux)*z1), 3)
            self.z[r[i]-1] = round((((1-cos0)*ux*uz + sin0*uy)*x1 + ((1-cos0)*uy*uz + sin0*ux)*y1 + (cos0 + (1-cos0)*(uz**2))*z1), 3)

    def rotate_list(self, r, c, start, b):
        if(c < len(r) and len(r) < self.topology.num_atoms):
            for i in range(start, len(self.topology.bonds), 2):
                if((r[c] == int(self.topology.bonds[i])) and (b != int(self.topology.bonds[i+1]))):
                    r.append(int(self.topology.bonds[i+1]))
                    start = i+2
                    break
                elif((r[c] == int(self.topology.bonds[i+1])) and (b != int(self.topology.bonds[i]))): 
                    r.append(int(self.topology.bonds[i]))
                    start = i+2
                    break
            for i in range(start, len(self.topology.bonds), 2):
                if((r[c] == int(self.topology.bonds[i])) and (b != int(self.topology.bonds[i+1]))):
                    r.append(int(self.topology.bonds[i+1]))
                elif((r[c] == int(self.topology.bonds[i+1])) and (b != int(self.topology.bonds[i]))): 
                    r.append(int(self.topology.bonds[i]))
            self.rotate_list(r, c+1, start, b)


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
           
           
molecule = Molecule()

pdb, psf, out_pdb, out_psf = get_cmd_line()

molecule.readPDB(pdb)
molecule.topology.readPSF(psf)
bond_list_C3 = molecule.bond_list(3)
angle_list_C3 = molecule.angle_list(3)
dihed_list_C3 = molecule.dihedral_list(3)
molecule.rotate(2, 3, np.pi/2)
molecule.writePDB("butane_opt_rotated.pdb")

for i in range(0, len(bond_list_C3), 2):
	print(bond_list_C3[i:i+2])
print("\n\n")
for i in range(0, len(angle_list_C3), 3):
	print(angle_list_C3[i:i+3])
print("\n\n")
for i in range(0, len(dihed_list_C3), 4):
	print(dihed_list_C3[i:i+4])
print("\n\n")
for i in range(0, len(molecule.topology.bonds), 2):
	print(molecule.topology.bonds[i:i+2])
