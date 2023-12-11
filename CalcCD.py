import numpy as np

from Functions import r_linear_interp
from numpy.linalg import norm
from numpy import dot


def calc_CD(ni: int, nj: int, vel_cd: object, CFL, Re, Pr, RES,
            cell_volume: object, cell_center: object,
            i_face_center: object, j_face_center: object,
            i_face_vector: object, j_face_vector: object):
    max_iter = 100000
    EPS = 1.0E-8
    ITER = 1
    T = np.zeros((ni + 1, nj + 1))
    T[0, :] = 1
    T[ni, :] = 2
    T_ex = T
    with open('Residuals.plt', 'w') as f:
        while ITER <= max_iter:
            RES = 0
            for i in range(1, ni):
                for j in range(1, nj):
                    dx_min = 100000
                    dx_max = 0
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

                        dnc = norm(cell_center[i_n, j_n, :] - cell_center[i, j, :])

                        nf = s_f / norm(s_f)
                        if 2 * dc < dx_min:
                            dx_min = 2 * dc
                        if 2 * dc > dx_max:
                            dx_min = 2 * dc

                        v_f = np.zeros(2)
                        v_f[0] = r_linear_interp(dc, dn, vel_cd[i, j, 0], vel_cd[i_n, j_n, 0])  # скорость на грани
                        v_f[1] = r_linear_interp(dc, dn, vel_cd[i, j, 1], vel_cd[i_n, j_n, 1])

                        tc = T_ex[i, j]
                        tn = T_ex[i_n, j_n]
                        tf = r_linear_interp(dc, dn, tc, tn)

                        if j_n == nj:
                            RES[i, j] = RES[i, j]
                        elif j_n == 0:
                            RES[i, j] += dot(nf, v_f) * tf * norm(s_f)
                        else:
                            RES[i, j] += (dot(nf, v_f) * tf - (tn - tc) / (dnc * Re * Pr)) * norm(s_f)

                    dx = 0.5 * (dx_min + dx_max)
                    dt = CFL * dx / 1
                    T[i, j] = T_ex[i, j] - dt * RES[i, j] / cell_volume[i-1, j-1]


            if ITER == 1:
                RESM = np.max(RES[1:ni-1, 1:nj-1])
            f.write(f'{ITER} {np.max(RES[1:ni-1,1:nj-1])/RESM}')

            ITER += 1
            T_ex = T

            if np.max(np.abs(RES[1:ni-1,1:nj-1])) < EPS:
                break

    return T


