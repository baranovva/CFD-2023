import Functions as funcs
import numpy as np

from OutputFields import output_fields
from CalcMetric import CalcMetric
from CalcDivergence import calc_divergence
from CalcGradient import calc_gradient
from CalcLaplacian import calc_laplacian

ni, nj, x, y = funcs.mesh_reader(file_name='base.msh')

# Calculate metric
cell_volume, cell_center, i_face_center, j_face_center, i_face_vector, j_face_vector = CalcMetric(ni, nj, x, y).run()

# Initiate field P
p = np.zeros((ni + 1, nj + 1))
grad_p_exact = np.zeros((ni + 1, nj + 1, 2))
for i in range(ni + 1):
    for j in range(nj + 1):
        p[i, j] = 5 * cell_center[i, j, 0] + 3 * cell_center[i, j, 1]
        grad_p_exact[:, :, 0] = 5
        grad_p_exact[:, :, 1] = 3

# funcs.figure(5, 'p', p[1:ni, 1:nj], 'turbo')

# Calculate gradient p
grad_p = np.zeros((ni + 1, nj + 1, 2))
n_iter = 10
for it in range(n_iter):
    grad_p = calc_gradient(ni=ni, nj=nj, p=p, grad_p=grad_p,
                           cell_volume=cell_volume,
                           cell_center=cell_center,
                           i_face_center=i_face_center,
                           j_face_center=j_face_center,
                           i_face_vector=i_face_vector,
                           j_face_vector=j_face_vector)

    grad_p_error = np.abs(1 - (grad_p / grad_p_exact))

print(f'Maximum GradPx-error: {np.max(grad_p_error[1:ni, 1:nj, 0])}')
print(f'Maximum GradPy-error: {np.max(grad_p_error[1:ni, 1:nj, 1])}')
# funcs.figure(5, 'grad_p_x', grad_p[1:ni, 1:nj, 0], 'turbo')
# funcs.figure(5, 'grad_p_y', grad_p[1:ni, 1:nj, 1], 'turbo')

# Initiate field V
v = np.zeros((ni + 1, nj + 1, 2))
div_v_p_exact = np.zeros((ni + 1, nj + 1))
for i in range(ni + 1):
    for j in range(nj + 1):
        v[i, j, 0] = 1 + cell_center[i, j, 0]
        v[i, j, 1] = 1 + cell_center[i, j, 1]
        div_v_p_exact[i, j] = 2 + 3 * cell_center[i, j, 0] + 3 * cell_center[i, j, 1]

# funcs.figure(5, 'vx', v[1:ni, 1:nj, 0], 'turbo')
# funcs.figure(5, 'vy', v[1:ni, 1:nj, 1], 'turbo')
funcs.figure(5, 'div_v_p_exact', div_v_p_exact[1:ni, 1:nj], 'turbo')

# Calculate divergence V and P
mode = 1
div_v_p = calc_divergence(mode=mode, ni=ni, nj=nj,
                          v=v, grad_p=grad_p, p=p,
                          cell_volume=cell_volume,
                          cell_center=cell_center,
                          i_face_center=i_face_center,
                          j_face_center=j_face_center,
                          i_face_vector=i_face_vector,
                          j_face_vector=j_face_vector)
div_v_p_error = np.abs(1 - (div_v_p / div_v_p_exact))

print(f'Maximum DivV-error: {np.max(div_v_p_error[1:ni, 1:nj])}')
funcs.figure(5, 'div_v', div_v_p[1:ni, 1:nj], 'turbo')
funcs.figure(5, 'div_v_error', div_v_p_error[1:ni, 1:nj], 'turbo')

# Calculate laplacian P
lap_p_exact = np.zeros((ni + 1, nj + 1))
for i in range(ni + 1):
    for j in range(nj + 1):
        lap_p_exact[i, j] = 4

lap_p = calc_laplacian(ni=ni, nj=nj, p=p, grad_p=grad_p,
                       cell_volume=cell_volume,
                       cell_center=cell_center,
                       i_face_center=i_face_center,
                       j_face_center=j_face_center,
                       i_face_vector=i_face_vector,
                       j_face_vector=j_face_vector)
lap_p_error = np.abs(1 - (lap_p / lap_p_exact))

print(f'Maximum LapP-error: {np.max(lap_p_error[1:ni, 1:nj])}')
funcs.figure(5, 'lap_p', lap_p[1:ni, 1:nj], 'turbo')

# Output fields to file
# output_file = 'data.dat'
# output_fields(file_name=output_file, ni=ni, nj=nj,
#               x=x, y=y, p=p, v=v,
#               grad_p=grad_p, grad_p_error=grad_p_error)
