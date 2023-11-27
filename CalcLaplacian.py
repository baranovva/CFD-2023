import numpy as np

from Functions import r_linear_interp
from numpy.linalg import norm


def calc_laplacian(ni: int, nj: int, p: object, grad_p: object,
                   cell_volume: object, cell_center: object,
                   i_face_center: object, j_face_center: object,
                   i_face_vector: object, j_face_vector: object) -> object:
    lap_p = np.zeros((ni + 1, nj + 1))
    for i in range(1, ni):
        for j in range(1, nj):
            for i_face in (1, 2, 3, 4):
                if i_face == 1:
                    i_n = i - 1
                    j_n = j
                    r_f = i_face_center[i - 1, j - 1, :]  # центр грани
                    s_f = - i_face_vector[i - 1, j - 1, :]
                elif i_face == 2:
                    i_n = i + 1
                    j_n = j
                    r_f = i_face_center[i, j - 1, :]
                    s_f = i_face_vector[i, j - 1, :]
                elif i_face == 3:
                    i_n = i
                    j_n = j - 1
                    r_f = j_face_center[i - 1, j - 1, :]
                    s_f = - j_face_vector[i - 1, j - 1, :]
                else:
                    i_n = i
                    j_n = j + 1
                    r_f = j_face_center[i - 1, j, :]
                    s_f = j_face_vector[i - 1, j, :]

                dc = norm(r_f[:] - cell_center[i, j, :])  # расстояние от границы до центра текущей ячейки,
                # радиус вектор центра ячейка, индексы пробегаем
                dn = norm(r_f[:] - cell_center[i_n, j_n, :])  # расстояние от границы до центра соседних ячеек
                # радиус вектор центра ячейка, индексы заграничные
                dnc = norm(cell_center[i_n, j_n, :] - cell_center[i, j, :])  # расстояние между центрами ячеек
                n_f = np.zeros_like(s_f)
                n_f[:] = s_f[:] / norm(s_f[:])  # eдиничный вектор внешней нормали

                # for skew correction
                rnc = (cell_center[i_n, j_n, :] - cell_center[i, j, :]) / dnc  # вектор по линии центров
                g_f = np.zeros_like(s_f)
                g_f[0] = r_linear_interp(dc, dn, grad_p[i, j, 0], grad_p[i_n, j_n, 0])
                g_f[1] = r_linear_interp(dc, dn, grad_p[i, j, 1], grad_p[i_n, j_n, 1])

                dpdn = (p[i_n, j_n] - p[i, j]) / dnc  # производная  dp/dn

                if dn < 1e-5:
                    dpdn_c = np.dot(grad_p[i, j, :], n_f[:])
                    dpdn = (5 * dpdn - 2 * dpdn_c) / 3  # 2 - ORDER
                    g_f = grad_p[i, j, :]

                dpdn += np.dot(n_f[:] - rnc[:], g_f[:])
                lap_p[i, j] += dpdn * norm(s_f[:])

            lap_p[i, j] /= cell_volume[i - 1, j - 1]
    return lap_p
