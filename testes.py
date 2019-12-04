import numpy as np
import sys
from itertools import chain
import matplotlib.pyplot as plt

#Liste todos os ângulos em que o atomo 3 é o centro.
#
#Resposta:
#
#    2 3 4
#    2 3 10
#    2 3 11

#Liste todos os diedros em que o atomo 3 está envolvido no centro da torção
#
#Resposta:
# 
#    1 2 3 4
#    1 2 3 10
#    1 2 3 11
#    
#    8 2 3 4
#    8 2 3 10
#    8 2 3 11
#    
#    9 2 3 4
#    9 2 3 10
#    9 2 3 11
#    
#    2 3 4 12
#    2 3 4 13
#    2 3 4 14
#    
#    10 3 4 12
#    10 3 4 13
#    10 3 4 14
#    
#    11 3 4 12
#    11 3 4 13
#    11 3 4 14

#Liste todos os atomos que devem ser rodados na torcao 2-3, do lado do atomo 3.
#
#Resposta: 3, 4, 10, 11, 12, 13, 14

def iter_data(fileobj):
    for line in fileobj:
        yield from line.split()

def read_triangular_array(path):
    with open(path) as fileobj:
        first = fileobj.readline().split()
        n = len(first)
        count = int(n*(n+1)/2)
        data = chain(first, iter_data(fileobj))
        return np.fromiter(data, int, count=count)

#def funcao1(tam):
#    l1 = np.fromiter(count = tam, dtype = int)
#    return l1
#    # ineficiente

def funcao2(tam):
    l2 = [None]*tam
    return l2

def funcao3(tam):
    l3 = np.zeros(tam, int)
    return l3
    
def funcao4(tam):
    l4 = np.empty(tam, int)
    return l4

def funcao5(tam):
    lista = [np.zeros(tam, int) for i in range(tam)]
    lista = [i for i in range(tam)]
    l5 = np.triu(lista)
    return lista, l5

# Para ver tempos: F10

tam = 4
#l1 = funcao1(tam*tam)
l2 = funcao2(tam*tam)
l3 = funcao3(tam*tam)
l4 = funcao4(tam*tam)
l, l5 = funcao5(tam)

d = dict()
X = ['c3', 'c3', 'c2', 'c2', 'c3']
Y = ['c3', 'hc', 'c2', 'hc', 'c2']
A = ['A', 'E', 'I', 'O', 'U', '']
B = ['B', 'C', 'D', 'F', 'G']

V = [28.51788850870531,
     28.529129963494892,
     28.55690923860382,
     28.514383625156793,
     28.526121507510585,
     28.480691477893675,
     28.44553584206169,
     28.453143714287318,
     28.45459964991983,
     28.481648866417075,
     28.498614662816415,
     28.489015127262178,
     28.506857590549732,
     28.530053976735353,
     28.54835606031772,
     28.56028719701376,
     28.505717953898458,
     28.4598439271044,
     28.486264955182175,
     28.5349073962644,
     28.558806223814077,
     28.555779046117276,
     28.567439526991997,
     28.618162593748334,
     28.58193549515719,
     28.56495089275606,
     28.530626307491115,
     28.555572184146698,
     28.537852363569545,
     28.588031903973274,
     28.614885141174568,
     28.577917882889523,
     28.540469290514704,
     28.51824391136676,
     28.539005222976954,
     28.5326315818481,
     28.551849388407454,
     28.518366357287054,
     28.497860272275844,
     28.495930441492824,
     28.500579684378827,
     28.5363266383939,
     28.522208643652654,
     28.5195104929224,
     28.52227097369648,
     28.531260171318127,
     28.540260925590225,
     28.53480123964264,
     28.53431651869327,
     28.54782686349725,
     28.515820786754745,
     28.581502957846013,
     28.60501136990519,
     28.573970842657307,
     28.593597586107208,
     28.57260851003717,
     28.609189086804964,
     28.596904660882,
     28.580696242891577,
     28.626270258877756,
     28.636262436985685,
     28.660601255141653,
     28.635513386435917,
     28.66537229721089,
     28.628707380805196,
     28.63356756134587,
     28.668411099955524,
     28.639728269052917,
     28.67813144568356,
     28.6465508254643,
     28.66078161156592,
     28.64762501248959,
     28.64762501248959]

for i in range(len(V)):
    plt.plot(i, V[i])
plt.xlabel('Iteracao')
plt.ylabel('Potencial')
plt.show()