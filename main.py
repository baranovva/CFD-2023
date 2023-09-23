from CalcMetric import CalcMetric
from time import perf_counter

import Functions
import CalcGradient
import OutputFields
import numpy as np


def mesh_reader(file_name):
    with open(file_name, 'r') as f:
        # Read nodes number (NI,NJ) from file with mesh
        ni, nj = map(int, f.readline().split())
        x = y = np.empty((ni, nj))
        # Read mesh from file
        for j in range(nj):
            for i in range(ni):
                x[i, j], y[i, j] = map(float, f.readline().split())

        return ni, nj, x, y


time_start = perf_counter()

mesh_file = 'base.msh'
NI, NJ, X, Y = mesh_reader(file_name=mesh_file)

# Calculate metric
CellVolume, CellCenter = CalcMetric(NI, NJ, X, Y).run()

# Initiate fields
P = np.empty((NI + 1, NJ + 1))
for j in range(NJ + 1):
    for i in range(NI + 1):
        P[i, j] = Functions.Pressure(CellCenter[i, j, 0], CellCenter[i, j, 1])

# Calculate gradient
grad_p = CalcGradient.B_CalcGradient(NI, NJ, P)

# Output fields to file
output_file = 'data.dat'
OutputFields.B_OutputFields(output_file, NI, NJ, X, Y, P, grad_p)

time_elapsed = perf_counter() - time_start
print(time_elapsed)
