'''
@author: Eleonora Pippia
'''

import numpy as np

class py_sapo:
    def __init__(self, nvar):
        self.nvar = nvar
        if self.nvar == 2:
            self.L = np.array([[1, 0], [0, 1], [1, 1], [1, -1]], dtype=np.double)
            self.T = np.array([[0, 1], [2, 3]], dtype=np.double)
            self.n_bundle = 2
            self.n_dir = 4
        elif self.nvar == 3:
            self.L = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1],
                               [1, 1, 0], [1, -1, 0], [0, 1, 1],
                               [0, 1, -1], [1, 0, 1], [1, 0, -1],
                               [1, 1, 1], [1, 1, -1], [1, -1, 1],
                               [1, -1, -1]], dtype=np.double)
            self.T = np.array([[0, 1, 2], [9, 3, 7], [10, 0, 4], [8, 11, 12], [5, 7, 6]], dtype=np.double)
            self.n_bundle = 5
            self.n_dir = 13
        self.offp = np.zeros(self.n_dir, dtype=np.double)
        self.offm = np.zeros(self.n_dir, dtype=np.double)

    def LTmatrix(self):
        return self.L, self.T, self.n_bundle

    def offset(self, cdd, ncons):
        if self.nvar == 2:
            for i in range(ncons):
                for j in range(self.n_dir):
                    if (self.L[j][0] == cdd[i][1]) & (self.L[j][1] == cdd[i][2]):
                        self.offm[j] = cdd[i][0]
                    elif (self.L[j][0] == -cdd[i][1]) & (self.L[j][1] == -cdd[i][2]):
                        self.offp[j] = cdd[i][0]
        elif self.nvar == 3:
            for i in range(ncons):
                for j in range(self.n_dir):
                    if (self.L[j][0] == cdd[i][1]) & \
                            (self.L[j][1] == cdd[i][2]) & \
                            (self.L[j][2] == cdd[i][3]):
                        self.offm[j] = cdd[i][0]
                    elif (self.L[j][0] == -cdd[i][1]) & \
                            (self.L[j][1] == -cdd[i][2]) & \
                            (self.L[j][2] == -cdd[i][3]):
                        self.offp[j] = cdd[i][0]

    def createRegions(self):
        regions = []
        self.outputbounds = np.zeros([pow(3, self.nvar), 2*self.n_dir], dtype=np.double)
        if self.nvar == 2:
            lb = np.array([-self.offm[0], -self.offm[1]], dtype=np.double)
            ub = np.array([self.offp[0], self.offp[1]], dtype=np.double)
            lbx1 = [[-lb[0], 1, 0], [2, 1, 0], [-2, 1, 0]]
            ubx1 = [[-2, -1, 0], [2, -1, 0], [ub[0], -1, 0]]
            lbx2 = [[-lb[1], 0, 1], [2, 0, 1], [-2, 0, 1]]
            ubx2 = [[-2, 0, -1], [2, 0, -1], [ub[1], 0, -1]]
            for i in range(3):
                for j in range(3):
                    regions.append(lbx1[j])
                    regions.append(ubx1[j])
                    regions.append(lbx2[2 - i])
                    regions.append(ubx2[2 - i])
                    self.outputbounds[j + i*3, :] = [1, 1, 2, 2, 1, 1, 2, 2]
                    if i == 0:
                        if j == 0:
                            self.outputbounds[j + i * 3, :] = np.concatenate((-self.L@[-1, 1], self.L@[-1, 1]), axis=0)
                        elif j == 2:
                            self.outputbounds[j + i * 3, :] = np.concatenate((-self.L @ [1, 1], self.L @ [1, 1]),
                                                                             axis=0)
                    elif i == 2:
                        if j == 0:
                            self.outputbounds[j + i * 3, :] = np.concatenate((-self.L @ [-1, -1], self.L @ [-1, -1]), axis=0)
                        elif j == 2:
                            self.outputbounds[j + i * 3, :] = np.concatenate((-self.L @ [1, -1], self.L @ [1, -1]),
                                                                             axis=0)

        elif self.nvar == 3:
            lb = np.array([-self.offm[0], -self.offm[1], -self.offm[2]], dtype=np.double)
            ub = np.array([self.offp[0], self.offp[1], self.offp[2]], dtype=np.double)
            lbx1 = [[-lb[0], 1, 0, 0], [2, 1, 0, 0], [-2, 1, 0, 0]]
            ubx1 = [[-2, -1, 0, 0], [2, -1, 0, 0], [ub[0], -1, 0, 0]]
            lbx2 = [[-lb[1], 0, 1, 0], [2, 0, 1, 0], [-2, 0, 1, 0]]
            ubx2 = [[-2, 0, -1, 0], [2, 0, -1, 0], [ub[1], 0, -1, 0]]
            lbx3 = [[-lb[2], 0, 0, 1], [2, 0, 0, 1], [-2, 0, 0, 1]]
            ubx3 = [[-2, 0, 0, -1], [2, 0, 0, -1], [ub[2], 0, 0, -1]]
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        regions.append(lbx1[k])
                        regions.append(ubx1[k])
                        regions.append(lbx2[j])
                        regions.append(ubx2[j])
                        regions.append(lbx3[i])
                        regions.append(ubx3[i])
                        self.outputbounds[k + j * 3 + i * 9, :] = [1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 1, 1, 1, 2, 2, 2, 2,
                                                           2, 2, 3, 3, 3, 3]
                        if i == 0:
                            if j == 0:
                                if k == 0:
                                    self.outputbounds[k + j * 3 + i * 9, :] = np.concatenate(
                                        (-self.L @ [-1, -1, -1], self.L @ [-1, -1, -1]), axis=0)
                                elif k == 2:
                                    self.outputbounds[k + j * 3 + i * 9, :] = np.concatenate(
                                        (-self.L @ [-1, -1, 1], self.L @ [-1, -1, 1]), axis=0)
                            elif j == 2:
                                if k == 0:
                                    self.outputbounds[k + j * 3 + i * 9, :] = np.concatenate(
                                        (-self.L @ [-1, 1, -1], self.L @ [-1, 1, -1]), axis=0)
                                elif k == 2:
                                    self.outputbounds[k + j * 3 + i * 9, :] = np.concatenate(
                                        (-self.L @ [-1, 1, 1], self.L @ [-1, 1, 1]), axis=0)
                        elif i == 2:
                            if j == 0:
                                if k == 0:
                                    self.outputbounds[k + j * 3 + i * 9, :] = np.concatenate(
                                        (-self.L @ [1, -1, -1], self.L @ [1, -1, -1]), axis=0)
                                elif k == 2:
                                    self.outputbounds[k + j * 3 + i * 9, :] = np.concatenate(
                                        (-self.L @ [1, -1, 1], self.L @ [1, -1, 1]), axis=0)
                            elif j == 2:
                                if k == 0:
                                    self.outputbounds[k + j * 3 + i * 9, :] = np.concatenate(
                                        (-self.L @ [1, 1, -1], self.L @ [1, 1, -1]), axis=0)
                                elif k == 2:
                                    self.outputbounds[k + j * 3 + i * 9, :] = np.concatenate(
                                        (-self.L @ [1, 1, 1], self.L @ [1, 1, 1]), axis=0)
        return regions

    def emptyoutputcons(self):
        output_cons = np.concatenate((self.L, -self.L), axis=0)
        output_cons = np.concatenate((np.zeros([self.n_dir*2, 1], dtype=np.double), output_cons), axis=1)
        return output_cons

    def comput_valOutputcons(self, cons):
        return self.outputbounds[cons, :]
