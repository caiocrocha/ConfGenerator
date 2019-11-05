class Molecule:
    def __init__(self):
        self.chain = list()
    class Chain:
            def __init__(self):
                self.id = list()    # "ATOM " or "HETATM"
                self.atom = list()  # atom name
                self.residue = list()   # residue name
                self.res_num = list()   # residue sequence number
                self.x = list() # orthogonal coordinates for X (in Angstroms)
                self.y = list() # orthogonal coordinates for Y (in Angstroms)
                self.z = list() # orthogonal coordinates for Z (in Angstroms)

    def readPDB(self, filename):
        with open(filename, "r") as inf:
            i = 0
            self.chain.append(self.Chain())
            for line in inf:
                line = line.split()
                if(line[0] == 'ATOM' or line[0] == 'HETATM'):
                    self.chain[i].id.append(line[0])
                    self.chain[i].atom.append(line[2])
                    self.chain[i].residue.append(line[3])
                    self.chain[i].res_num.append(int(line[4]))
                    self.chain[i].x.append(float(line[5]))
                    self.chain[i].y.append(float(line[6]))
                    self.chain[i].z.append(float(line[7]))
                elif(line[0] == 'TER'):
                    i += 1
                    self.chain.append(self.Chain())
                    
    def writePDB(self, filename):
        with open(filename, "w+") as outf:
            for i in range(len(self.chain)):
                for j in range(len(self.chain[i].id)):
                    outf.write("{:6s}{:5d} {:^5s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n".format(
                    self.chain[i].id[j],
                    j+1,    # atom serial number
                    self.chain[i].atom[j],
                    self.chain[i].residue[j],
                    '', # alternate location indicator
                    self.chain[i].res_num[j],
                    '', # chain identifier
                    self.chain[i].x[j],
                    self.chain[i].y[j],
                    self.chain[i].z[j],
                    1.00, 0.00, '', ''))    # occupancy, temperature factor, element symbol, charge on the atom
                if(j < len(self.chain[i].id)):
                    outf.write("TER\n")
            outf.write("END\n")
    
molecule1 = Molecule()
#molecule2 = Molecule()
   
molecule1.readPDB("butane.pdb")
#molecule2.readPDB("ethane.pdb")

molecule1.writePDB("butane1.pdb")
#molecule2.writePDB("ethane1.pdb")

