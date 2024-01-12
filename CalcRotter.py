import numpy as np

from Functions import interp
from numpy.linalg import norm
from numba import njit, prange

@njit(cache=True, parallel=True)
def calc_rotter(ni: int, nj: int, v: object,
                cell_volume: object, cell_center: object,
                i_face_center: object, j_face_center: object,
                i_face_vector: object, j_face_vector: object):
    rot_v = np.zeros((ni + 1, nj + 1))
    for i in prange(1, ni):
        for j in prange(1, nj):
            for i_face in prange(1, 5):
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

                dc = norm(r_f[:] - cell_center[i, j, :])  # расстояние от границы до центра текущей ячейки,
                # радиус вектор центра ячейка, индексы пробегаем
                dn = norm(r_f[:] - cell_center[i_n, j_n, :])  # расстояние от границы до центра соседних ячеек
                # радиус вектор центра ячейка, индексы заграничные

                v_f = np.empty(2)
                v_f[0] = interp(dc, dn, v[i, j, 0], v[i_n, j_n, 0])  # скорость на грани
                v_f[1] = interp(dc, dn, v[i, j, 1], v[i_n, j_n, 1])

                rot_v[i, j] += s_f[0] * v_f[1] - s_f[1] * v_f[0]

            rot_v[i, j] /= cell_volume[i, j]

    return rot_v
