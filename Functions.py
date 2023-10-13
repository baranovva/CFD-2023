import numpy as np


def pressure(x: object, y: object) -> object:
    return 5 * x + 3 * y


def r_linear_interp(d1, d2, x1, x2) -> float:
    return (x1 * d2 + x2 * d1) / (d1 + d2)


def grad_p_exact(
        ni: int,
        nj: int,
        x: object,
        y: object
) -> object:
    grad_p_ex = np.empty((ni + 1, nj + 1, 2))
    grad_p_ex[:, :, 0] = 5
    grad_p_ex[:, :, 1] = 3

    return grad_p_ex


def velocity(
        ni: int,
        nj: int,
        x: object,
        y: object
) -> object:
    v = np.empty((ni + 1, nj + 1, 2))
    v[:, :, 0] = 1 + x
    v[:, :, 1] = 1 + y

    return v


def div_velocity_exact(x: object, y: object):
    return 1


def div_velocity_p_exact(x: object, y: object):
    return 2 + 3 * x + 3 * y


def calc_rot_exact(x: object, y: object):
    return 2


def r_laplacian_p_exact(x: object, y: object):
    return 4
