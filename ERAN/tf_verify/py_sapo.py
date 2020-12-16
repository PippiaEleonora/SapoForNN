'''
@author: Eleonora Pippia
'''

import numpy as np
import math

class py_sapo:
    def __init__(self, nvar, n_dir, cdd):
        self.nvar = nvar
        self.n_dir = n_dir
        self.L = np.zeros([n_dir, nvar], dtype=np.double)
        self.offp = np.zeros(self.n_dir, dtype=np.double)
        self.offm = np.zeros(self.n_dir, dtype=np.double)
        relation = -np.ones(2*n_dir, dtype=np.int)
        for i in range(n_dir * 2 - 1):
            if relation[i] == -1:
                for j in range(i + 1, n_dir * 2):
                    if np.prod(cdd[j, 1:] == -cdd[i, 1:]):
                        relation[i] = j
                        relation[j] = i
        count = 0
        for i in range(n_dir):
            while relation[count] < count:
                count = count + 1

            self.L[i, :] = cdd[count, 1:]
            self.offm[i] = cdd[count][0]
            self.offp[i] = cdd[relation[count]][0]
            count = count + 1

        directions = list(range(n_dir))
        selected = []
        index = np.zeros(nvar, dtype=np.int)
        T = []
        self.n_bundle = 0
        lineardependent = 0

        while (len(directions) > nvar-1) and (not lineardependent):
            for i in range(nvar):
                index[i] = directions[i]
            count = nvar-1
            matrix = self.L[index, :]
            while np.linalg.det(matrix) == 0:
                count = count + 1
                if count < len(directions):
                    index[nvar-1] = directions[count]
                    matrix = self.L[index, :]
                else:
                    lineardependent = 1
                    break

            if not lineardependent:
                T.append(index.tolist())
                for i in range(nvar):
                    directions.remove(index[i])
                    selected.append(index[i])
                self.n_bundle = self.n_bundle + 1

        while len(directions) > 0:
            for i in range(min(len(directions), nvar-1)):
                index[i] = directions[i]
                id1 = i
            id1 = id1 + 1
            for i in range(id1):
                directions.remove(index[i])
            for j in range(id1, nvar):
                index[j] = selected[j-len(directions)]
            id2 = nvar-len(directions)
            count = id2-1
            matrix = self.L[[index], :]
            while np.linalg.det(matrix) == 0:
                count = count + 1
                if count < len(selected):
                    index[nvar-1] = selected[count]
                    matrix = self.L[index, :]
                else:
                    print('Error!\n')

            T.append(index.tolist())
            self.n_bundle = self.n_bundle + 1
        self.T = np.asarray(T, dtype=np.double)


    def LTmatrix(self):
        return self.L, self.T, self.n_bundle

    def createRegions(self):
        regions = []
        self.outputbounds = np.zeros([pow(3, self.nvar), 2*self.n_dir], dtype=np.double)
        if self.nvar == 2:
            lbpoint = np.array([-1, -1, 1], dtype=np.double)
            ubpoint = np.array([-1, 1, 1], dtype=np.double)
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
                    regions.append(lbx2[i])
                    regions.append(ubx2[i])

                    points = [[x, y] for x in [lbpoint[j], ubpoint[j]] for y in [lbpoint[i], ubpoint[i]]]
                    points = np.asarray(points, dtype=np.double)
                    points = np.unique(points, axis=0)

                    result1 = np.min(points @ np.transpose(self.L), axis=0)
                    result2 = np.max(points @ np.transpose(self.L), axis=0)
                    self.outputbounds[j + i * 3, :] = self.outputbounds[j + i * 3, :] = np.concatenate((-result1, result2), axis=0)

        elif self.nvar == 3:
            lbpoint = np.array([-1, -1, 1], dtype=np.double)
            ubpoint = np.array([-1, 1, 1], dtype=np.double)
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

                        points = [[x, y, z] for x in [lbpoint[k], ubpoint[k]] for y in [lbpoint[j], ubpoint[j]]
                                  for z in [lbpoint[i], ubpoint[i]]]
                        points = np.asarray(points, dtype=np.double)
                        points = np.unique(points, axis=0)

                        result1 = np.min(points@np.transpose(self.L), axis=0)
                        result2 = np.max(points @ np.transpose(self.L), axis=0)
                        self.outputbounds[k + j * 3 + i * 9, :] = self.outputbounds[k + j * 3 + i * 9, :] = np.concatenate(
                            (-result1, result2), axis=0)
        return regions

    def emptyoutputcons(self):
        output_cons = np.concatenate((self.L, -self.L), axis=0)
        output_cons = np.concatenate((np.zeros([self.n_dir*2, 1], dtype=np.double), output_cons), axis=1)
        return output_cons

    def comput_valOutputcons(self, cons):
        return self.outputbounds[cons, :]
