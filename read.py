class Molecula:
    def __init__(self):
        self.atomo = list()
        self.residuo = list()
        self.num_residuo = list()
        self.x = list()
        self.y = list()
        self.z = list()

molecula = Molecula()

with open("butane.pdb", "r") as inf:
    for line in inf:
        line = line.split()
        if(line[0] == 'ATOM'):
            molecula.atomo.append(line[2])
            molecula.residuo.append(line[3])
            molecula.num_residuo.append(line[4])
            molecula.x.append(line[5])
            molecula.y.append(line[6])
            molecula.z.append(line[7])

with open("saida.pdb", "w+") as outf:
    for i in range(len(molecula.atomo)):
        outf.write("ATOM\t")
        outf.write(str(i))
        outf.write('\t')
        outf.write(str(molecula.atomo[i]))
        outf.write('\t')
        outf.write(str(molecula.residuo[i]))
        outf.write('\t')
        outf.write(str(molecula.num_residuo[i]))
        outf.write('\t')
        outf.write(str(molecula.x[i]))
        outf.write('\t')
        outf.write(str(molecula.y[i]))
        outf.write('\t')
        outf.write(str(molecula.z[i]))
        outf.write("\t1.00\t0.00\n")
    outf.write("TER\nEND\n")