import numpy as np
import matplotlib.pyplot as plt


def mesh_reader(file_name: str) -> object:
    with (open(file_name, 'r') as f):
        # Read nodes number (ni,nj) from file with mesh
        ni, nj = map(int, f.readline().split())
        x = np.empty((ni, nj))
        y = np.empty((ni, nj))
        # Read mesh from file
        for i in range(ni):
            for j in range(nj):
                x[i, j], y[i, j] = map(float, f.readline().split())

        return ni, nj, x, y


def figure(figsize, title, data, cmap):
    plt.figure(figsize=(figsize, figsize))
    plt.imshow(data, cmap=cmap)
    plt.xlabel("I")
    plt.ylabel("J")
    plt.title(title)
    plt.show()


def r_linear_interp(d1, d2, x1, x2) -> float:
    return (x1 * d2 + x2 * d1) / (d1 + d2)
