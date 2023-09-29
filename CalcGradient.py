import numpy as np

from numpy.linalg import norm


def calc_gradient(
        ni, nj, p,
        cell_volume,
        cell_center,
        i_face_center,
        j_face_center,
        i_face_vector,
        j_face_vector
):
    grad_p = np.empty((ni + 1, nj + 1, 2))
    r_c = r_n = np.zeros(2)

    for i in range(ni - 1):
        for j in range(nj - 1):
            gp = np.zeros(2)

            for i_face in range(1, 5):

                if i_face == 1:
                    i_n = i - 1
                    j_n = j
                    r_f = i_face_center[i, j, :]
                    s_f = - i_face_vector[i, j, :]

                elif i_face == 2:
                    i_n = i + 1
                    j_n = j
                    r_f = i_face_center[i + 1, j, :]
                    s_f = i_face_vector[i + 1, j, :]

                elif i_face == 3:
                    i_n = i
                    j_n = j - 1
                    r_f = j_face_center[i, j, :]
                    s_f = - j_face_vector[i, j, :]

                else:
                    i_n = i
                    j_n = j + 1
                    r_f = j_face_center[i, j + 1, :]
                    s_f = j_face_vector[i, j + 1, :]

                vol = cell_volume[i, j]
                r_c[:] = cell_center[i, j, :]
                r_n[:] = cell_center[i_n, j_n, :]

                dc = norm(r_f[:] - r_c[:])
                dn = norm(r_f[:] - r_n[:])

            grad_p[i, j, :] = gp[:] / vol
    return grad_p
