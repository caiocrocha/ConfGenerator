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

molecule1 = Topology()
#molecule2 = Topology()

molecule1.readPSF("butane.psf")
#molecule2.readPSF("ethane.psf")

molecule1.writePSF("butane1.psf")
#molecule2.writePDB("ethane1.psf")
