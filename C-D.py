import Functions as funcs
import numpy as np

from CalcCD import calc_cd
from OutputFields import output_fields

ni, nj, x, y = funcs.mesh_reader(file_name='InputMesh.txt')

with (open('VelocityField.txt', 'r') as f):
    _ = f.readline()
    v_cd = np.array([[[float(x) for x in f.readline().split()] for _ in range(ni + 1)] for _ in range(nj + 1)])

with (open('Temperature.txt', 'r') as f):
    _ = f.readline()
    T_e = np.array([[float(f.readline()) for _ in range(ni + 1)] for _ in range(nj + 1)])

T = calc_cd(ni=ni, nj=nj, v_cd=v_cd, x=x, y=y, max_iter=10, CFL=0.1, Re=100, Pr=1, first_scheme=False)

T_error = np.abs(1 - (T / T_e))
print(f'Maximum T error: {np.max(T_error[1:ni, 1:nj])}')
funcs.figure(5, 'T', T[1:ni, 1:nj], 'turbo')
funcs.figure(5, 'T_error', T_error[1:ni, 1:nj], 'turbo')

output_fields(file_name='data_cd.dat', ni=ni, nj=nj, x=x, y=y, v=v_cd, t=T, t_error=T_error)
