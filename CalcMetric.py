import numpy as np
from numpy.linalg import norm


class CalcMetric:
    def __init__(self, ni: int, nj: int, x: object, y: object):
        self.NI = ni
        self.NJ = nj
        self.X = x
        self.Y = y
        self.i_face_center = np.empty((self.NI + 1, self.NJ, 2))
        self.i_face_vector = np.empty((self.NI + 1, self.NJ, 2))
        self.j_face_center = np.empty((self.NI, self.NJ + 1, 2))
        self.j_face_vector = np.empty((self.NI, self.NJ + 1, 2))

    def face_centers_vectors_i(self):
        for i in range(1, self.NI + 1):
            for j in range(1, self.NJ):
                self.i_face_vector[i, j, 0] = self.Y[i, j + 1] - self.Y[i, j]
                self.i_face_vector[i, j, 1] = -(self.X[i, j + 1] - self.X[i, j])
                self.i_face_center[i, j, 0] = 0.5 * (self.X[i, j] + self.X[i, j + 1])
                self.i_face_center[i, j, 1] = 0.5 * (self.Y[i, j] + self.Y[i, j + 1])

    def face_centers_vectors_j(self):
        for i in range(1, self.NI):
            for j in range(1, self.NJ + 1):
                self.j_face_vector[i, j, 0] = - (self.Y[i + 1, j] - self.Y[i, j])
                self.j_face_vector[i, j, 1] = self.X[i + 1, j] - self.X[i, j]
                self.j_face_center[i, j, 0] = 0.5 * (self.X[i, j] + self.X[i + 1, j])
                self.j_face_center[i, j, 1] = 0.5 * (self.Y[i, j] + self.Y[i + 1, j])

    def cell_volumes(self) -> object:
        cell_volume = np.empty((self.NI, self.NJ))
        r = np.empty(2)
        for i in range(1, self.NI):
            for j in range(1, self.NJ):
                r[0] = self.X[i + 1, j + 1] - self.X[i, j]
                r[1] = self.Y[i + 1, j + 1] - self.Y[i, j]
                cell_volume[i, j] = 0.5 * (np.dot(self.i_face_vector[i, j, :], r[:]) +
                                           np.dot(self.j_face_vector[i, j, :], r[:]))

        return cell_volume

    def cell_centers(self) -> object:
        cell_center = np.empty((self.NI + 1, self.NJ + 1, 2))

        # FOR INNER CELLS: CENTER OF CONTOUR
        for i in range(1, self.NI):
            for j in range(1, self.NJ):
                numerator = (self.i_face_center[i, j, :] * norm(self.i_face_vector[i, j, :]) +
                             self.i_face_center[i + 1, j, :] * norm(self.i_face_vector[i + 1, j, :]) +
                             self.j_face_center[i, j, :] * norm(self.j_face_vector[i, j, :]) +
                             self.j_face_center[i, j + 1, :] * norm(self.j_face_vector[i, j + 1, :]))

                denominator = (norm(self.i_face_vector[i, j, :]) + norm(self.i_face_vector[i + 1, j, :]) +
                               norm(self.j_face_vector[i, j, :]) + norm(self.j_face_vector[i, j + 1, :]))

                cell_center[i, j, :] = numerator / denominator

        # FOR DUMMY CELLS ON BOUNDARIES: CELL CENTER = FACE CENTER
        for j in range(1, self.NJ):  # I-BOUNDARIES
            cell_center[0, j, :] = self.i_face_center[1, j, :]
            cell_center[self.NI, j, :] = self.i_face_center[self.NI, j, :]

        for i in range(1, self.NI):  # J-BOUNDARIES
            cell_center[i, 0, :] = self.j_face_center[i, 1, :]
            cell_center[i, self.NJ, :] = self.j_face_center[i, self.NJ, :]

        return cell_center

    def run(self) -> object:
        self.face_centers_vectors_i()
        self.face_centers_vectors_j()
        cell_center = self.cell_centers()
        cell_volume = self.cell_volumes()

        return cell_volume, cell_center, self.i_face_center, self.j_face_center, self.i_face_vector, self.j_face_vector


def read_metrics(file_name: str, ni: int, nj: int) -> object:
    cell_center = np.empty((ni + 1, nj + 1, 2))
    cell_volume = np.empty((ni, nj))
    i_face_center = np.empty((ni + 1, nj, 2))
    i_face_vector = np.empty((ni + 1, nj, 2))
    j_face_center = np.empty((ni, nj + 1, 2))
    j_face_vector = np.empty((ni, nj + 1, 2))
    with (open(file_name, 'r') as f):
        for i in range(ni + 1):
            for j in range(nj + 1):
                cell_center[i, j, 0], cell_center[i, j, 1] = map(float, f.readline().split())
        for i in range(1, ni):
            for j in range(1, nj):
                cell_volume[i, j] = float(f.readline())
        for i in range(1, ni + 1):
            for j in range(1, nj):
                i_face_center[i, j, 0], i_face_center[i, j, 1] = map(float, f.readline().split())
        for i in range(1, ni + 1):
            for j in range(1, nj):
                i_face_vector[i, j, 0], i_face_vector[i, j, 1] = map(float, f.readline().split())
        for i in range(1, ni):
            for j in range(1, nj + 1):
                j_face_center[i, j, 0], j_face_center[i, j, 1] = map(float, f.readline().split())
        for i in range(1, ni):
            for j in range(1, nj + 1):
                j_face_vector[i, j, 0], j_face_vector[i, j, 1] = map(float, f.readline().split())

        return cell_volume, cell_center, i_face_center, j_face_center, i_face_vector, j_face_vector
