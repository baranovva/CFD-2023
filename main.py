import Functions as funcs
import numpy as np
import time

from CalcMetric import CalcMetric
from CalcGradient import calc_gradient
from CalcDivergence import calc_divergence
from CalcLaplacian import calc_laplacian
from CalcRotter import calc_rotter
from OutputFields import output_fields
from joblib import Parallel, delayed

ni, nj, x, y = funcs.mesh_reader(file_name='input/base.msh')  # base.msh fine.msh skew.msh
cell_volume, cell_center, i_face_center, j_face_center, i_face_vector, j_face_vector = CalcMetric(ni, nj, y, x).run()


# Initiate field P and calculate gradient p
def calculate_p(i, j, cell_center):
    p_ij = cell_center[i, j, 0] ** 2 + cell_center[i, j, 1] ** 2 + cell_center[i, j, 0] + cell_center[i, j, 1]
    grad_p_exact_ij = np.array([2 * cell_center[i, j, 0] + 1, 2 * cell_center[i, j, 1] + 1])
    return (i, j, p_ij, grad_p_exact_ij)


p = np.zeros((ni + 1, nj + 1))
grad_p_exact = np.zeros((ni + 1, nj + 1, 2))

start_time = time.time()
results = Parallel(n_jobs=1)(delayed(calculate_p)(i, j, cell_center) for i in range(ni + 1) for j in range(nj + 1))
print(f'Time:{time.time() - start_time} s')

for result in results:
    i, j, p_ij, grad_p_exact_ij = result
    p[i, j] = p_ij
    grad_p_exact[i, j, :] = grad_p_exact_ij

grad_p = np.zeros((ni + 1, nj + 1, 2))
n_iter = 10
start_time = time.time()
for _ in range(n_iter):
    grad_p = calc_gradient(ni=ni, nj=nj, p=p, grad_p=grad_p,
                           cell_volume=cell_volume, cell_center=cell_center,
                           i_face_center=i_face_center, j_face_center=j_face_center,
                           i_face_vector=i_face_vector, j_face_vector=j_face_vector)
print(f'Time:{time.time() - start_time} s')
grad_p_error = np.abs((grad_p_exact - grad_p) / grad_p_exact)
print(f'Maximum GradPx-error: {np.max(grad_p_error[1:ni, 1:nj, 0])}')
print(f'Maximum GradPy-error: {np.max(grad_p_error[1:ni, 1:nj, 1])}')

# Initiate field V and calculate divergence V and P
v = np.empty((ni + 1, nj + 1, 2))
div_v_p_exact = np.empty((ni + 1, nj + 1))
v[:, :, 0] = 1 + cell_center[:, :, 0]
v[:, :, 1] = 1 + cell_center[:, :, 1]
div_v_p_exact[:, :] = (4 * (cell_center[:, :, 0] ** 2 + cell_center[:, :, 1] ** 2) +
                       5 * (cell_center[:, :, 0] + cell_center[:, :, 1]) + 2)

mode = 1
div_v_p = calc_divergence(mode=mode, ni=ni, nj=nj, v=v, grad_p=grad_p, p=p,
                          cell_volume=cell_volume, cell_center=cell_center,
                          i_face_center=i_face_center, j_face_center=j_face_center,
                          i_face_vector=i_face_vector, j_face_vector=j_face_vector)
div_v_p_error = np.abs(1 - (div_v_p / div_v_p_exact))
print(f'Maximum DivV_P-error: {np.max(div_v_p_error[1:ni, 1:nj])}')

# Calculate laplacian P
lap_p_exact = np.zeros((ni + 1, nj + 1))
lap_p_exact[:, :] = 4

lap_p = calc_laplacian(ni=ni, nj=nj, p=p, grad_p=grad_p,
                       cell_volume=cell_volume, cell_center=cell_center,
                       i_face_center=i_face_center, j_face_center=j_face_center,
                       i_face_vector=i_face_vector, j_face_vector=j_face_vector)
lap_p_error = np.abs(1 - (lap_p / lap_p_exact))
print(f'Maximum LapP-error: {np.max(lap_p_error[1:ni, 1:nj])}')

# Calculate rotV
rot_v_exact = np.empty((ni + 1, nj + 1))
rot_v_exact[:, :] = 0

rot_v = calc_rotter(ni=ni, nj=nj, v=v,
                    cell_volume=cell_volume, cell_center=cell_center,
                    i_face_center=i_face_center, j_face_center=j_face_center,
                    i_face_vector=i_face_vector, j_face_vector=j_face_vector)
rot_v_error = np.abs(1 - (rot_v / rot_v_exact))
print(f'Maximum rotV error: {np.max(rot_v_error[1:ni, 1:nj])}')

# Output fields to file
output_file = 'output/data.dat'
output_fields(file_name=output_file, ni=ni, nj=nj, x=x, y=y, p=p, v=v,
              grad_p=grad_p, grad_p_error=grad_p_error, div_v_p=div_v_p, div_v_p_error=div_v_p_error)
