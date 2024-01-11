import numpy as np

from CalcMetric import CalcMetric
from Functions import r_linear_interp
from numpy.linalg import norm


def calc_cd(ni: int, nj: int, v_cd: object, max_iter: int, x: object, y: object,
            CFL=0.001, Re=100, Pr=1, first_scheme=True):
    (cell_volume, cell_center, i_face_center,
     j_face_center, i_face_vector, j_face_vector) = CalcMetric(ni, nj, y, x).run()

    EPS = 1e-7
    T = np.ones((ni + 1, nj + 1))
    T[:, nj] = 2
    T_ex = T
    res = np.array([])
    v_f = np.zeros(2)

    for iteration in range(1, max_iter + 1):
        RES = np.zeros((ni + 1, nj + 1))
        for i in range(1, ni):
            for j in range(1, nj):
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

                    dc = norm(r_f[:] - cell_center[i, j, :])  # расстояние от границы до центра текущей ячейки,
                    # радиус вектор центра ячейка, индексы пробегаем
                    dn = norm(r_f[:] - cell_center[i_n, j_n, :])  # расстояние от границы до центра соседних ячеек
                    # радиус вектор центра ячейка, индексы заграничные

                    v_f[0] = r_linear_interp(dc, dn, v_cd[i, j, 0], v_cd[i_n, j_n, 0])  # скорость на грани
                    v_f[1] = r_linear_interp(dc, dn, v_cd[i, j, 1], v_cd[i_n, j_n, 1])
                    tc = T_ex[i, j]
                    tn = T_ex[i_n, j_n]

                    if first_scheme:  # central
                        tf = r_linear_interp(dc, dn, tc, tn)
                    else:
                        if np.dot(s_f, v_f) >= 0:
                            tf = tc
                        elif dn < 1e-6:
                            tf = 2 * tn - tc
                        else:
                            tf = tn

                    nf = s_f / norm(s_f)

                    if j_n == 0:
                        RES[i, j] += np.dot(nf, v_f) * tf * norm(s_f)
                    elif j_n != nj:
                        dnc = norm(cell_center[i_n, j_n, :] - cell_center[i, j, :])
                        RES[i, j] += (np.dot(nf, v_f) * tf - (tn - tc) / (dnc * Re * Pr)) * norm(s_f)

                T[i, j] = T_ex[i, j] - CFL * RES[i, j] / np.sqrt(cell_volume[i, j])
                if T[i, j] < 1:
                    T[i, j] = 1

        if iteration == 1: res_max = np.max(RES[1:ni - 1, 1:nj - 1])
        res = np.append(res, [iteration, np.max(RES[1:ni - 1, 1:nj - 1]) / res_max])
        if np.max(np.abs(RES[1:ni - 1, 1:nj - 1])) < EPS: break
        iteration += 1
        T_ex = T

    with open('Residuals.plt', 'w') as f:
        np.savetxt(f, res.reshape(-1, 2), fmt='%25.17f')
    return T
