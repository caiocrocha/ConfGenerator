filename='butane.psf'
bond_list=[]
angle_list=[]
dihedral_list=[]

f = open(filename,"r")
line = f.readline()

while line:

	if ("!NBOND" in line) :
		nbonds=int(line.split()[0])
		nlines=int(nbonds/4)

		if(nbonds%2 != 0):
			nlines+=1

		for x in range(nlines):
			bond_list.extend(f.readline().split())

	if ("!NTHETA" in line) :
		nangles=int(line.split()[0])
		nlines=int(nangles/3)

		if(nangles%2 != 0):
			nlines+=1

		for x in range(nlines):
			angle_list.extend(f.readline().split())

	if ("!NPHI" in line) :
		ndihedrals=int(line.split()[0])
		nlines=int(ndihedrals/2)

		if(ndihedrals%2 != 0):
			nlines+=1

		for x in range(nlines):
			dihedral_list.extend(f.readline().split())


	line = f.readline()		

# Tarefas:
# Liste todos os atomos ligados ao atomo 3.
atom_list=[]
for i in range(0,len(bond_list),2):
	if (int(bond_list[i])==3):
		atom_list.append(int(bond_list[i+1]))
	elif (int(bond_list[i+1])==3): 
		atom_list.append(int(bond_list[i]))

print(atom_list)
'''
# Liste todos os ângulos em que o atomo 3 é o centro.
2 3 4
2 3 10
2 3 11
# Liste todos os diedros em que o atomo 3 está envolvido no centro da torção
resposta 
1 2 3 4
1 2 3 10
1 2 3 11
8 2 3 4
8 2 3 10
8 2 3 11
9 2 3 4
9 2 3 10
9 2 3 11
2 3 4 12
2 3 4 13
2 3 4 14
10 3 4 12
10 3 4 13
10 3 4 14
11 3 4 12
11 3 4 13
11 3 4 14


# Liste todos os atomos que devem ser rodados na torcao 2-3, do lado do atomo 3.
resposta 3,4,10,11,12,13,14
'''
