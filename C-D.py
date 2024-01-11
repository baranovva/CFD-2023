import Functions as funcs
import numpy as np

from CalcMetric import CalcMetric
from CalcGradient import calc_gradient
from CalcDivergence import calc_divergence
from CalcLaplacian import calc_laplacian
from CalcRotter import calc_rotter
from OutputFields import output_fields
from CalcCD import calc_CD

# Calculate C-D
nk, nl, x_cd, y_cd = funcs.mesh_reader(file_name='InputMesh.txt')
(cell_volume_cd, cell_center_cd, i_face_center_cd,
 j_face_center_cd, i_face_vector_cd, j_face_vector_cd) = CalcMetric(nk, nl, x_cd, y_cd).run()


with (open('VelocityField.txt', 'r') as f):
    nk_temp, nl_temp = map(int, f.readline().split())
    v_cd = np.zeros((nk_temp, nl_temp, 2))
    for i in range(nk_temp):
        for j in range(nl_temp):
            v_cd[i, j, 0], v_cd[i, j, 1] = map(float, f.readline().split())

with (open('Temperature.txt', 'r') as f):
    nk_temp, nl_temp = map(int, f.readline().split())
    T_e = np.zeros((nk_temp, nl_temp))
    for i in range(nk_temp):
        for j in range(nl_temp):
            T_e[i, j] = float(f.readline())

# funcs.figure(5, 'v_cd', v_cd[:, :, 1], 'turbo')
# funcs.figure(5, 'T', T[:, :], 'turbo')

CFL = 0.01
Re = 100
Pr = 1
T = calc_CD(ni=nk, nj=nl, vel_cd=v_cd,
            CFL=CFL, Re=Re, Pr=Pr,
            cell_volume=cell_volume_cd,
            cell_center=cell_center_cd,
            i_face_center=i_face_center_cd,
            j_face_center=j_face_center_cd,
            i_face_vector=i_face_vector_cd,
            j_face_vector=j_face_vector_cd)

T_error = np.abs(1 - (T / T_e))
print(f'Maximum T error: {np.max(T_error[:, :])}')
# funcs.figure(5, 'T', T[1:ni, 1:nj], 'turbo')