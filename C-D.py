import Functions as funcs
import numpy as np
import time

from CalcCD import calc_cd
from OutputFields import output_fields

ni, nj, x, y = funcs.mesh_reader(file_name='input/InputMesh.txt')  # InputMesh.txt Mesh.msh

with (open('input/VelocityField.txt', 'r') as f):  # VelocityField.txt Vel_Central.txt Vel_FOU.txt
    _ = f.readline()
    v_cd = np.array([[[float(x) for x in f.readline().split()] for _ in range(nj + 1)] for _ in range(ni + 1)])

with (open('input/Temperature.txt', 'r') as f):  # Temperature.txt T_Central.txt T_FOU.txt
    _ = f.readline()
    T_e = np.array([[float(f.readline()) for _ in range(nj + 1)] for _ in range(ni + 1)])

t = time.time()
T = calc_cd(ni=ni, nj=nj, v_cd=v_cd, x=x, y=y, max_iter=10, CFL=0.01, Re=100, central_scheme=True)
print(time.time() - t)

T_error = np.abs(1 - (T / T_e))
print(f'Maximum T error: {np.max(T_error[1:ni, 1:nj])}')
output_fields(file_name='cd/data_cd.dat', ni=ni, nj=nj, x=x, y=y, v=v_cd, t=T, t_error=T_error)
funcs.figure(5, 'T', T[:, :], 'turbo')
funcs.figure(5, 'T_error', T_error[:, :], 'turbo')
