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
            flag = ''
            for line in inf:
                if(line != '\n'):
                    line = line.split()
                    if('!' in line[1]):
                        flag = line[1]
                    if('!NATOM' in line[1]):
                        self.num_atoms = int(line[0])
                    elif('!NBOND' in line[1]):
                        self.num_bonds = int(line[0])
                    elif('!NTHETA' in line[1]):
                        self.num_angles = int(line[0])
                    elif('!NPHI' in line[1]):
                        self.num_dihedrals = int(line[0])
                    elif('!NATOM' in flag):
                        self.segid.append(line[1])
                        self.resid.append(int(line[2]))
                        self.resname.append(line[3])
                        self.name.append(line[4])
                        self.type.append(line[5])
                        self.charge.append(float(line[6]))
                        self.mass.append(float(line[7]))
                    elif('!NBOND' in flag):
                        self.bonds.extend(list(zip(*[iter(line)] * 2)))
                    elif('!NTHETA' in flag):
                        self.angles.extend(list(zip(*[iter(line)] * 3)))
                    elif('!NPHI' in flag):
                        self.dihedrals.extend(list(zip(*[iter(line)] * 4)))

    def writePSF(self, filename):
        with open(filename, "w+") as outf:
            outf.write("PSF CHEQ EXT XPLOR\n\n{:10d} !NTITLE\n\n\n{:10d} !NATOM\n".format(1, self.num_atoms))
            for i in range(self.num_atoms):
                outf.write("{:10d} {:5s}{:5d} {:>10s}{:^10s}  {:5s}{:14.6f} {:14.4f}\n".format(
                i+1,
                self.segid[i],
                self.resid[i],
                self.resname[i],
                self.name[i],
                self.type[i],
                self.charge[i],
                self.mass[i]))
            outf.write("{:10d} !NBOND\n".format(self.num_bonds))
            for i in range(self.num_bonds):
                outf.write("{:>10s}{:>10s}".format(self.bonds[i][0], self.bonds[i][1]))
                if(i%4 == 0 and i != 0):
                    outf.write("\n")
            outf.write("{:10d} !NTHETA\n".format(self.num_angles))
            for i in range(self.num_angles):
                outf.write("{:>10s}{:>10s}".format(self.angles[i][0], self.angles[i][1]))
                if(i%3 == 0 and i != 0):
                    outf.write("\n")
            outf.write("{:10d} !NBOND\n".format(self.num_dihedrals))
            for i in range(self.num_dihedrals):
                outf.write("{:>10s}{:>10s}".format(self.dihedrals[i][0], self.dihedrals[i][1]))
                if(i%2 == 0 and i != 0):
                    outf.write("\n")

molecule1 = Topology()
#molecule2 = Topology()

molecule1.readPSF("butane.psf")
#molecule2.readPSF("ethane.psf")

molecule1.writePSF("butane1.psf")
#molecule2.writePDB("ethane1.psf")
