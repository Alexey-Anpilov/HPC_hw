import numpy as np
import random
import struct
n = 1000

M = [random.random() for i in range(n * n)]

test_m1 = [1, -1, 1,
           1, 0, 1,
           1, 1, 2]

test_m2 = [1, 0, 0,
           0, 0, 1,
           0, 1, 0]

matrices = {"matrix.dat" : M, "test_matrix1.dat" : test_m1, "test_matrix2.dat": test_m2}

for m_name in matrices.keys():
    with open(m_name, 'wb') as file:
        b = struct.pack('%sd' % len(matrices[m_name]), *matrices[m_name])
        file.write(b)
