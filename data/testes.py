# Remark line goes here
# MASS
# c3 12.010        0.878
# hc 1.008         0.135

# BOND
# c3-c3  300.90 K  1.538 b0
# c3-hc  330.60   1.097

# ANGLE
# c3-c3-c3   62.860 K angulo    111.510 angulo eq (graus)
# c3-c3-hc   46.340     109.800
# hc-c3-hc   39.400     107.580

# DIHE
# c3-c3-c3-c3   1    0.180 K        0.000 angulo         -3.000 multiplicidade
# c3-c3-c3-c3   1    0.250       180.000          -2.000
# c3-c3-c3-c3   1    0.200       180.000           1.000
# c3-c3-c3-hc   1    0.160         0.000           3.000
# hc-c3-c3-hc   1    0.150         0.000           3.000

# IMPROPER

# NONBON
#   c3          1.9080  0.1094
#   hc          1.4870  0.0157

import sys
import numpy as np
from time import time
import argparse
def get_cmd_line():
    parser = argparse.ArgumentParser(description='profiler')
    parser.add_argument('--tam', 	  action='store', dest='tam',	   required=True, help='tamanho do vetor')
    arguments = parser.parse_args()
    # Assign from cmd line.
    return arguments.tam
tam = int(get_cmd_line())
lista0 = np.array([np.zeros(tam, dtype=bool), np.zeros(tam, dtype=bool)])
lista1 = np.zeros(tam, dtype=bool)
lista2 = np.zeros(tam, dtype='uint8')
lista3 = []
for i in range(1, tam, 3):
	lista0[0][i] = True
	lista0[0][i-1] = True
	lista1[i] = True
	lista2[i] = 1
	lista1[i-1] = True
	lista2[i-1] = 1
t0 = time()
for i in range(1, tam):
	if lista0[0][i] and not lista0[1][i] and lista0[0][i] != lista0[0][i-1]:
		lista0[1][i] = True
t00 = time() - t0
t1 = time()
for i in range(1, tam):
	if lista1[i] and lista1[i] != lista1[i-1] and i not in lista3:
		lista3.append(i)
t01 = time() - t1
t2 = time()
for i in range(1, tam):
	if lista2[i] and lista1[i] != lista1[i-1] and lista2[i] == 1:
		lista2[i] = 2
t02 = time() - t2
print('{:.16f} {:.16f} {:.16f}'.format(float(t00), float(t01), float(t02)))
