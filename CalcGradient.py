import numpy as np


def B_CalcGradient(ni, nj, p):
    grad_p = np.zeros((ni + 1, nj + 1, 2))
    for i in range(ni - 1):
        for j in range(nj - 1):
            grad_p[i, j, :] = p[i, j]
    return grad_p
