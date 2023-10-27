import numpy as np

from Functions import r_linear_interp
from numpy.linalg import norm
from numpy import dot


def calc_divergence(ni: int, nj: int, v: object,
                    cell_volume: object,
                    cell_center: object,
                    i_face_center: object,
                    j_face_center: object,
                    i_face_vector: object,
                    j_face_vector: object):
    div_v = np.zeros((ni + 1, nj + 1))
    r_c = np.zeros(2)
    r_n = np.zeros(2)

    for i in range(2, ni - 2):
        for j in range(2, nj - 2):
            for i_face in (1, 2, 3, 4):
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

                r_c[:] = cell_center[i, j, :]  # радиус вектор центра ячейка, индексы пробегаем
                r_n[:] = cell_center[i_n, j_n, :]  # радиус вектор центра ячейка, индексы заграничные

                dc = norm(r_f[:] - r_c[:])  # расстояние от границы до центра текущей ячейки
                dn = norm(r_f[:] - r_n[:])  # расстояние от границы до центра соседних ячеек

                v_f = np.zeros(2)
                v_f[0] = r_linear_interp(dc, dn, v[i, j, :], v[i_n, j_n, :])  # скорость на грани
                v_f[1] = r_linear_interp(dc, dn, v[i, j, :], v[i_n, j_n, :])

                div_v[i, j] += dot(v_f[:] * s_f[:])

            div_v[i, j] = div_v[:] / cell_volume[i, j]

    return div_v
