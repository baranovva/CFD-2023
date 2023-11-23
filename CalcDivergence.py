import numpy as np

from Functions import r_linear_interp
from numpy.linalg import norm
from numpy import dot


def calc_divergence(mode: int, ni: int, nj: int, v: object, grad_p: object, p: object,
                    cell_volume: object, cell_center: object,
                    i_face_center: object, j_face_center: object,
                    i_face_vector: object, j_face_vector: object):
    div_v = np.zeros((ni + 1, nj + 1))
    for i in range(1, ni):
        for j in range(1, nj):
            for i_face in range(1, 5):
                if i_face == 1:
                    i_n = i - 1
                    j_n = j
                    r_f = i_face_center[i - 1, j - 1, :]
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

                v_f = np.zeros(2)
                v_f[0] = r_linear_interp(dc, dn, v[i, j, 0], v[i_n, j_n, 0])  # скорость на грани
                v_f[1] = r_linear_interp(dc, dn, v[i, j, 1], v[i_n, j_n, 1])

                if mode == 1:
                    p_f = r_linear_interp(dc, dn, p[i, j], p[i_n, j_n])
                elif mode == 2:  # противопоточная 1 порядка
                    if dot(s_f, v_f) >= 0:
                        p_f = p[i, j]
                    elif dn <= 1e-6:  # экстраполяция на грани
                        p_f = 2 * p[i_n, j_n] - p[i, j]
                    else:
                        p_f = p[i_n, j_n]
                else:  # противопоточная 2 порядка
                    if dot(s_f, v_f) >= 0:
                        p_f = p[i, j] + dot(r_f[:] - cell_center[i, j, :], grad_p[i, j, :])
                    elif dn <= 1e-6:  # экстраполяция на грани
                        pn = 2 * p[i_n, j_n] - p[i, j]
                        gc = dot(grad_p[i, j, :], cell_center[i, j, :] - r_f[:])
                        gb = p[i, j] - p[i_n, j_n]
                        gn = 4 * gb - 3 * gc
                        p_f = pn + gn
                    else:
                        p_f = p[i_n, j_n] + dot(r_f[:] - cell_center[i_n, j_n, :], grad_p[i_n, j_n, :])

                div_v[i, j] += p_f * dot(v_f[:], s_f[:])

            div_v[i, j] = div_v[i, j] / cell_volume[i - 1, j - 1]

    return div_v
