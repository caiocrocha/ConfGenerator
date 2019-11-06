import numpy as np

class Molecule:
    def __init__(self):
        self.topology = self.Topology()
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

    def rotate(self, a, b, axis, angle):   # bond a-b, rotate chain from b // angle in radians
        r = list()
        for i in range(self.topology.num_bonds):
            if(a not in self.topology.bonds[i]):
                if(b == int(self.topology.bonds[i][0])):
                    r.append(int(self.topology.bonds[i][1]))
                elif(b == int(self.topology.bonds[i][1])):
                    r.append(int(self.topology.bonds[i][0]))
        for j in range(len(r)):
            for i in range(self.topology.num_bonds):
                if(b not in self.topology.bonds[i]):
                    if(r[j] == int(self.topology.bonds[i][0])):
                        r.append(int(self.topology.bonds[i][1]))
                    elif(r[j] == int(self.topology.bonds[i][1])):
                        r.append(int(self.topology.bonds[i][0]))
        for i in range(len(r)):
            print(r[i])
            '''
        x1 = self.x[b]
        y1 = self.y[b]
        z1 = self.z[b]
        if(axis == 'x' or axis == 'X'):
            self.y[b] = y1*np.cos(angle) - z1*np.sin(angle)
            self.z[b] = y1*np.sin(angle) + z1*np.cos(angle)
        elif(axis == 'y' or axis == 'Y'):
            self.x[b] = x1*np.cos(angle) - z1*np.sin(angle)
            self.z[b] = x1*np.sin(angle) + z1*np.cos(angle)
        elif(axis == 'z' or axis == 'Z'):
            self.y[b] = x1*np.cos(angle) - y1*np.sin(angle)
            self.z[b] = x1*np.sin(angle) + y1*np.cos(angle)
        for i in range(len(r)):
            x1 = self.x[r[i]]
            y1 = self.y[r[i]]
            z1 = self.z[r[i]]
            if(axis == 'x' or axis == 'X'):
                self.y[r[i]] = y1*np.cos(angle) - z1*np.sin(angle)
                self.z[r[i]] = y1*np.sin(angle) + z1*np.cos(angle)
            elif(axis == 'y' or axis == 'Y'):
                self.x[r[i]] = x1*np.cos(angle) - z1*np.sin(angle)
                self.z[r[i]] = x1*np.sin(angle) + z1*np.cos(angle)
            elif(axis == 'z' or axis == 'Z'):
                self.y[r[i]] = x1*np.cos(angle) - y1*np.sin(angle)
                self.z[r[i]] = x1*np.sin(angle) + y1*np.cos(angle)
                '''

molecule = Molecule()

molecule.readPDB("butane.pdb")
molecule.topology.readPSF("butane.psf")

molecule.rotate(2, 3, 'x', np.pi/2)
