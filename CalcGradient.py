import numpy as np

from Functions import interp
from numpy.linalg import norm
from numba import njit, prange

@njit(cache=True, parallel=True)
def calc_gradient(ni: int, nj: int, p: object, grad_p: object,
                  cell_volume: object, cell_center: object,
                  i_face_center: object, j_face_center: object,
                  i_face_vector: object, j_face_vector: object) -> object:
    for i in prange(1, ni):
        for j in prange(1, nj):
            gp = np.zeros(2)
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

                p_e = interp(dc, dn, p[i, j], p[i_n, j_n])  # давление на грани

                # коордиаты до точки e
                r_e = np.zeros(2)
                r_e[0] = interp(dc, dn, cell_center[i, j, 0], cell_center[i_n, j_n, 0])
                r_e[1] = interp(dc, dn, cell_center[i, j, 1], cell_center[i_n, j_n, 1])

                # grad в точке е
                gp_e = np.zeros(2)
                gp_e[0] = interp(dc, dn, grad_p[i, j, 0], grad_p[i_n, j_n, 0])
                gp_e[1] = interp(dc, dn, grad_p[i, j, 1], grad_p[i_n, j_n, 1])

                # лин интерполяция радиус вектора в точке e
                p_f = p_e + np.dot(r_f[:] - r_e[:], gp_e[:])

                gp[:] += p_f * s_f[:]  # gradp

            grad_p[i, j, :] = gp[:] / cell_volume[i, j]
    return grad_p
