import numpy as np

from CalcMetric import CalcMetric
from Functions import interp
from numpy.linalg import norm
from numpy import dot


def calc_cd(ni: int, nj: int, v_cd: object, max_iter: int, x: object,
            y: object, CFL=0.01, Re=100, Pr=1, central_scheme=True):
    (cell_volume, cell_center, i_face_center,
     j_face_center, i_face_vector, j_face_vector) = CalcMetric(ni, nj, y, x).run()

    T = np.ones((ni + 1, nj + 1))
    T[:, nj] = 2  # T[:, nj]
    T_ex = np.ones_like(T)
    res = np.array([])
    v_f = np.empty(2)
    cfl_sqrt_volume_v = CFL / np.sqrt(cell_volume[1:ni, 1:nj] * norm(v_cd[1:ni, 1:nj, :], axis=2))
    Re_Pr = Re * Pr

    for iteration in range(max_iter):
        T_ex[:, :] = T[:, :]
        RES = np.zeros_like(T)
        for i in range(1, ni):
            for j in range(1, nj):
                for face in range(1, 5):
                    if face == 1:
                        i_n = i - 1
                        j_n = j
                        r_f = i_face_center[i, j, :]
                        s_f = - i_face_vector[i, j, :]
                    elif face == 2:
                        i_n = i + 1
                        j_n = j
                        r_f = i_face_center[i + 1, j, :]
                        s_f = i_face_vector[i + 1, j, :]
                    elif face == 3:
                        i_n = i
                        j_n = j - 1
                        r_f = j_face_center[i, j, :]
                        s_f = - j_face_vector[i, j, :]
                    else:
                        i_n = i
                        j_n = j + 1
                        r_f = j_face_center[i, j + 1, :]
                        s_f = j_face_vector[i, j + 1, :]
                    if not (i_n == ni or i_n == 0):  # i_n == ni or i_n == 0
                        dc = norm(r_f[:] - cell_center[i, j, :])  # расстояние от границы до
                        # центра текущей ячейки, вектор центра ячейка, индексы пробегаем
                        dn = norm(r_f[:] - cell_center[i_n, j_n, :])  # расстояние от границы до
                        # центра соседних ячеек радиус, вектор центра ячейка, индексы заграничные

                        v_f[0] = interp(dc, dn, v_cd[i, j, 1], v_cd[i_n, j_n, 1])  # скорость на грани
                        v_f[1] = interp(dc, dn, v_cd[i, j, 0], v_cd[i_n, j_n, 0])

                        if central_scheme:  # central
                            tf = interp(dc, dn, T_ex[i, j], T_ex[i_n, j_n])
                        else:  # fou
                            if dot(s_f[:], v_f[:]) >= 0:
                                tf = T_ex[i, j]
                            elif dn < 1e-6:
                                tf = 2 * T_ex[i_n, j_n] - T_ex[i, j]
                            else:
                                tf = T_ex[i_n, j_n]

                        RES[i, j] += ((dot(s_f[:] / norm(s_f[:]), v_f[:])) * tf - (T_ex[i_n, j_n] - T_ex[i, j])
                                      / (norm(cell_center[i_n, j_n, :] - cell_center[i, j, :])
                                         * Re_Pr)) * norm(s_f[:])

        T[1:ni, 1:nj] = T_ex[1:ni, 1:nj] - RES[1:ni, 1:nj] * cfl_sqrt_volume_v[:, :]

        max_res_for_iteration = np.abs(np.max((RES[1:ni, 1:nj])))
        res = np.append(res, [iteration, max_res_for_iteration])
        if max_res_for_iteration < 1e-14:
            break

        iteration += 1

    max_res = res[1]
    res = res.reshape(-1, 2)
    res[:, 1] = res[:, 1] / max_res
    with open('output/Residuals.plt', 'w') as f:
        np.savetxt(f, res, fmt='%25.17f')

    return T
