# butane_opt.pdb
# CRYST1    0.000    0.000    0.000  90.00  90.00  90.00 P 1           1
# ATOM      1  C1  BUT S   1       2.913   1.158  -1.990  0.00  0.00      SYS   
# ATOM      2  C2  BUT S   1       1.285   0.913  -1.606  0.00  0.00      SYS   
# ATOM      3  C3  BUT S   1       1.069  -0.008  -0.123  0.00  0.00      SYS   
# ATOM      4  C4  BUT S   1      -0.495  -0.311   0.434  0.00  0.00      SYS   
# ATOM      5  H1  BUT S   1       3.527   1.241  -1.075  0.00  0.00      SYS   
# ATOM      6  H2  BUT S   1       3.317   0.314  -2.570  0.00  0.00      SYS   
# ATOM      7  H3  BUT S   1       3.058   2.075  -2.573  0.00  0.00      SYS   
# ATOM      8  H4  BUT S   1       0.774   1.885  -1.517  0.00  0.00      SYS   
# ATOM      9  H5  BUT S   1       0.791   0.385  -2.436  0.00  0.00      SYS   
# ATOM     10  H6  BUT S   1       1.591  -0.972  -0.233  0.00  0.00      SYS   
# ATOM     11  H7  BUT S   1       1.616   0.513   0.679  0.00  0.00      SYS   
# ATOM     12  H8  BUT S   1      -0.542  -0.263   1.531  0.00  0.00      SYS   
# ATOM     13  H9  BUT S   1      -1.207   0.432   0.030  0.00  0.00      SYS   
# ATOM     14  H10 BUT S   1      -0.842  -1.300   0.113  0.00  0.00      SYS   
# END

# Liste todos os ângulos em que o atomo 3 é o centro.

# Resposta:

#     2 3 4
#     2 3 10
#     2 3 11

# Liste todos os diedros em que o atomo 3 está envolvido no centro da torção

# Resposta:

#     1 2 3 4
#     1 2 3 10
#     1 2 3 11
    
#     8 2 3 4
#     8 2 3 10
#     8 2 3 11
    
#     9 2 3 4
#     9 2 3 10
#     9 2 3 11
    
#     2 3 4 12
#     2 3 4 13
#     2 3 4 14
    
#     10 3 4 12
#     10 3 4 13
#     10 3 4 14
    
#     11 3 4 12
#     11 3 4 13
#     11 3 4 14

# Liste todos os atomos que devem ser rodados na torcao 2-3, do lado do atomo 3.

# Resposta: 3, 4, 10, 11, 12, 13, 14

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
