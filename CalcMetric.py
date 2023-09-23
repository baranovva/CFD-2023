import numpy as np

from numpy.linalg import norm


class CalcMetric:
    def __init__(self, ni, nj, x, y):
        self.NI = ni
        self.NJ = nj
        self.X = x
        self.Y = y
        self.r = np.empty(2)
        self.i_face_center = self.i_face_vector = np.empty((self.NI, self.NJ - 1, 2))
        self.j_face_center = self.j_face_vector = np.empty((self.NI - 1, self.NJ, 2))

    def face_centers_vectors_i(self):
        for i in range(self.NI):
            for j in range(self.NJ - 1):
                self.r[0] = self.X[i, j + 1] - self.X[i, j]
                self.r[1] = self.Y[i, j + 1] - self.Y[i, j]
                self.i_face_vector[i, j, 0] = self.r[1]
                self.i_face_vector[i, j, 1] = -self.r[0]
                self.i_face_center[i, j, 0] = 0.5 * (self.X[i, j] + self.X[i, j + 1])
                self.i_face_center[i, j, 1] = 0.5 * (self.Y[i, j] + self.Y[i, j + 1])

    def face_centers_vectors_j(self):
        for i in range(self.NI - 1):
            for j in range(self.NJ):
                self.r[0] = self.X[i + 1, j] - self.X[i, j]
                self.r[1] = self.Y[i + 1, j] - self.Y[i, j]
                self.j_face_vector[i, j, 0] = -self.r[1]
                self.j_face_vector[i, j, 1] = self.r[0]
                self.j_face_center[i, j, 0] = 0.5 * (self.X[i, j] + self.X[i + 1, j])
                self.j_face_center[i, j, 1] = 0.5 * (self.Y[i, j] + self.Y[i + 1, j])

    def cell_volumes(self):
        cell_volume = np.empty((self.NI - 1, self.NJ - 1))

        for i in range(self.NI - 1):
            for j in range(self.NJ - 1):
                self.r[0] = self.X[i + 1, j + 1] - self.X[i, j]
                self.r[1] = self.Y[i + 1, j + 1] - self.Y[i, j]
                cell_volume[i, j] = (0.5 * np.dot(self.i_face_vector[i, j], self.r) +
                                     0.5 * np.dot(self.j_face_vector[i, j], self.r))

        return cell_volume

    def cell_centers(self):
        cell_center = np.empty((self.NI + 1, self.NJ + 1, 2))

        # FOR INNER CELLS: CENTER OF CONTOUR
        for i in range(self.NI - 1):
            for j in range(self.NJ - 1):
                numerator = (self.i_face_center[i, j] * norm(self.i_face_vector[i, j]) +
                             self.i_face_center[i + 1, j] * norm(self.i_face_vector[i + 1, j]) +
                             self.j_face_center[i, j] * norm(self.j_face_vector[i, j]) +
                             self.j_face_center[i, j + 1] * norm(self.j_face_vector[i, j + 1]))

                denominator = (norm(self.i_face_vector[i, j]) + norm(self.i_face_vector[i + 1, j]) +
                               norm(self.j_face_vector[i, j]) + norm(self.j_face_vector[i, j + 1]))

                cell_center[i, j, :] = numerator / denominator

        # FOR DUMMY CELLS ON BOUNDARIES: CELL CENTER = FACE CENTER
        for n_bound in range(2):
            if n_bound == 0:
                i_bound = j_bound = 1
                i_out = j_out = 0
            else:
                i_bound = i_out = -1
                j_bound = j_out = -1

            for j in range(self.NJ - 2):  # I-BOUNDARIES
                cell_center[i_out, j, :] = self.i_face_center[i_bound, j, :]

            for i in range(self.NI - 2):  # J-BOUNDARIES
                cell_center[i, j_out, :] = self.j_face_center[i, j_bound, :]

        return cell_center

    def run(self):
        self.face_centers_vectors_i()
        self.face_centers_vectors_j()
        cell_center = self.cell_centers()
        cell_volume = self.cell_volumes()

        return cell_volume, cell_center
