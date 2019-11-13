class Molecule:
    def __init__(self):
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
    
molecule1 = Molecule()
#molecule2 = Molecule()
   
molecule1.readPDB("butane_opt.pdb")
#molecule2.readPDB("ethane.pdb")

molecule1.writePDB("butane_opt1.pdb")
#molecule2.writePDB("ethane1.pdb")