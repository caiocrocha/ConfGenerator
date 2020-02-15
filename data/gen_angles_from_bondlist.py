# bond_list=[1, 2, 1, 3, 2, 3]
bond_list = [1, 2, 1, 6, 2, 3, 3, 4, 4, 5, 5, 6]

# metodo burro
for i in range(0, len(bond_list), 2):
    bond_list.append(bond_list[i + 1])
    bond_list.append(bond_list[i])

angle_list = []

for i in range(0, len(bond_list), 2):

    for j in range(i + 2, len(bond_list), 2):

        if bond_list[j] == bond_list[i]:
            angle_list.append(bond_list[i + 1])
            angle_list.append(bond_list[i])
            angle_list.append(bond_list[j + 1])

#            if bond_list[j] == bond_list[i+1]:
#                angle_list.append(bond_list[i])
#                angle_list.append(bond_list[j])
#                angle_list.append(bond_list[j+1])

print('BONDS')
for i in range(0, int(len(bond_list) / 2), 2):
    print(bond_list[i], bond_list[i + 1])

print('ANGLES')
for i in range(0, len(angle_list), 3):
    print(angle_list[i], angle_list[i + 1], angle_list[i + 2])
