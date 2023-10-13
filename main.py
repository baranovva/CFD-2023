import Functions as F
import numpy as np

from CalcGradient import calc_gradient
from OutputFields import output_fields
from CalcMetric import CalcMetric


def mesh_reader(file_name: str) -> object:
    with open(file_name, 'r') as f:
        # Read nodes number (ni,nj) from file with mesh
        ni, nj = map(int, f.readline().split())
        x = y = np.empty((ni, nj))

        # Read mesh from file
        for i in range(ni):
            for j in range(nj):
                x[i, j], y[i, j] = map(float, f.readline().split())

        return ni, nj, x, y


mesh_file = 'base.msh'
ni, nj, x, y = mesh_reader(file_name=mesh_file)

# Calculate metric
(cell_volume,
 cell_center,
 i_face_center,
 j_face_center,
 i_face_vector,
 j_face_vector) = CalcMetric(ni, nj, x, y).run()

# Initiate fields
p = np.empty((ni + 1, nj + 1))
div_v_exact = rot_v_exact = lap_p_exact = np.empty((ni + 1, nj + 1))

for i in range(ni + 1):
    for j in range(nj + 1):
        p[i, j] = F.pressure(
                x=cell_center[i, j, 0],
                y=cell_center[i, j, 1]
        )

        grad_p_exact = F.grad_p_exact(
                ni, nj,
                x=cell_center[i, j, 0],
                y=cell_center[i, j, 1],
        )

        # v = F.velocity(
        #         ni, nj,
        #         x=cell_center[i, j, 1],
        #         y=cell_center[i, j, 2]
        # )
        #
        # div_v_exact[i, j] = F.div_velocity_exact(
        #         x=cell_center[i, j, 1],
        #         y=cell_center[i, j, 2]
        # )
        #
        # rot_v_exact[i, j] = F.calc_rot_exact(
        #         x=cell_center[i, j, 1],
        #         y=cell_center[i, j, 2]
        # )
        #
        # lap_p_exact[i, j] = F.r_laplacian_p_exact(
        #         x=cell_center[i, j, 1],
        #         y=cell_center[i, j, 2]
        # )

# Calculate gradient
n_iter = 10
for it in range(0, n_iter):
    grad_p = calc_gradient(
            ni=ni, nj=nj,
            p=p,
            cell_volume=cell_volume,
            cell_center=cell_center,
            i_face_center=i_face_center,
            j_face_center=j_face_center,
            i_face_vector=i_face_vector,
            j_face_vector=j_face_vector
    )
    grad_p_error = abs(1 - (grad_p / grad_p_exact))
    print('Maximum GradPx-error:', np.max(grad_p_error[1:ni, 1:nj, 0]))
    print('Maximum GradPy-error:', np.max(grad_p_error[1:ni, 1:nj, 1]))

# Calculate gradient
# div_v =

# div_v_error = abs(div_v_exact - div_v) / div_v_exact

# print('Maximum DivV-error:', np.max(div_v_error))


# Calculate rotor
# rot_v =

# rot_v_error = abs(rot_v_exact-rot_v) / rot_v_exact

# print('Maximum RotV-error:', np.max(rot_v_error))


# Calculate laplasian
# lap_p =

# lap_p_error = abs(lap_p_exact - lap_p) / lap_p_exact

# print('Maximum LapP-error:', np.max(lap_p_error))


# Output fields to file
output_file = 'data.dat'
output_fields(
        file_name=output_file,
        ni=ni,
        nj=nj,
        x=x,
        y=y,
        p=p,
        grad_p=grad_p,
        grad_p_error=grad_p_error
)
