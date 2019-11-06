class Molecule:
    def __init__(self):
        self.id = list()    # "ATOM " or "HETATM" (or "TER")
        self.atom = list()  # atom name
        self.residue = list()   # residue name
        self.res_num = list()   # residue sequence number
        self.x = list() # orthogonal coordinates for X (in Angstroms)
        self.y = list() # orthogonal coordinates for Y (in Angstroms)
        self.z = list() # orthogonal coordinates for Z (in Angstroms)

    def readPDB(self, filename):
        with open(filename, "r") as inf:
            for line in inf:
                line = line.split()
                if(line[0] == 'ATOM' or line[0] == 'HETATM'):
                    self.id.append(line[0])
                    self.atom.append(line[2])
                    self.residue.append(line[3])
                    self.res_num.append(int(line[4]))
                    self.x.append(float(line[5]))
                    self.y.append(float(line[6]))
                    self.z.append(float(line[7]))
                elif(line[0] == 'TER'):
                    self.id.append(line[0])
                    self.atom.append('')
                    self.residue.append('')
                    self.res_num.append('')
                    self.x.append('')
                    self.y.append('')
                    self.z.append('')
                    
    def writePDB(self, filename):
        with open(filename, "w+") as outf:
            for i in range(len(self.id)):
                if(self.id[i] == 'ATOM' or self.id[i] == 'HETATM'):
                    outf.write("{:6s}{:5d} {:^5s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(
                    self.id[i],
                    i+1,    # atom serial number
                    self.atom[i],
                    self.residue[i],
                    '', # alternate location indicator
                    self.res_num[i],
                    '', # chain identifier
                    self.x[i],
                    self.y[i],
                    self.z[i],
                    1.00, 0.00, '', ''))    # occupancy, temperature factor, element symbol, charge on the atom
                elif(self.id[i] == 'TER'):
                    outf.write("TER\n")
            outf.write("END\n")
    
molecule1 = Molecule()
#molecule2 = Molecule()
   
molecule1.readPDB("butane.pdb")
#molecule2.readPDB("ethane.pdb")

molecule1.writePDB("butane1.pdb")
#molecule2.writePDB("ethane1.pdb")

