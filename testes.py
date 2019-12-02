import numpy as np
import sys
from itertools import chain

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
