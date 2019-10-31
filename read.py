class Molecula:
    def __init__(self):
        self.atomo = list()
        self.residuo = list()
        self.num_residuo = list()
        self.x = list()
        self.y = list()
        self.z = list()

molecula = Molecula()

def read_pdb(filename)
	with open(filename, "r") as inf:
		for line in inf:
			line = line.split()
			if(line[0] == 'ATOM'):
				molecula.atomo.append(line[2])
				molecula.residuo.append(line[3])
				molecula.num_residuo.append(line[4])
				molecula.x.append(float(line[5]))
				molecula.y.append(float(line[6]))
				molecula.z.append(float(line[7]))

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

with open("saida2.pdb", "w+") as outf:
    for i in range(len(molecula.atomo)):
        outf.write("{:6s}{:5d} {:^4s}{:1s}{:3s} {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format('ATOM',
          i,
          molecula.atomo[i],
          molecula.residuo[i],
          molecula.num_residuo[i],
          '',
          molecula.x[i],
          molecula.y[i],
          molecula.z[i],
          1.00,0.00
          )
         )
    outf.write("TER\nEND\n")

