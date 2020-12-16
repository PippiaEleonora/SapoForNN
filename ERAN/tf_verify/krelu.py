from elina_scalar import *
from elina_dimension import *
from elina_linexpr0 import *
from elina_abstract0 import *
from fppoly import *

import numpy as np
import cdd
import time
import itertools
import multiprocessing
import math

from config import config
"""
From reference manual:
http://web.mit.edu/sage/export/cddlib-094b.dfsg/doc/cddlibman.ps

CDD H-representaion format:
each row represents b + Ax >= 0
example: 2*x_1 - 3*x_2 >= 1 translates to [-1, 2, -3]

CDD V-representaion format:
<code> x_1 x_2 ... x_d
code is 1 for extremal point, 0 for generating direction (rays)
example: extreme point (2, -1, 3) translates to [1, 2, -1, 3]
all polytopes generated here should be closed, hence code=1
"""
import pandas as pd
from ctypes import cdll
from ctypes import *
from py_sapo import py_sapo
import numpy.ctypeslib
from poly_approx import poly_approx
boolean_flag=True



def generate_linexpr0(offset, varids, coeffs):
    # returns ELINA expression, equivalent to sum_i(varids[i]*coeffs[i])
    assert len(varids) == len(coeffs)
    n = len(varids)

    linexpr0 = elina_linexpr0_alloc(ElinaLinexprDiscr.ELINA_LINEXPR_SPARSE, n)
    cst = pointer(linexpr0.contents.cst)
    elina_scalar_set_double(cst.contents.val.scalar, 0)

    for i, (x, coeffx) in enumerate(zip(varids, coeffs)):
        linterm = pointer(linexpr0.contents.p.linterm[i])
        linterm.contents.dim = ElinaDim(offset + x)
        coeff = pointer(linterm.contents.coeff)
        elina_scalar_set_double(coeff.contents.val.scalar, coeffx)

    return linexpr0


class Krelu:
    def __init__(self, cdd_hrepr):
        array_2d_double = np.ctypeslib.ndpointer(dtype=np.uintp, ndim=1, flags='C')
        #nikos: poly approximation
        global boolean_flag
        if config.poly_dynamic is False:
            sapolib = cdll.LoadLibrary("../../Sapo/libsapo_dyn_lib.so")
            sapolib.computeSapo_small.argtypes = [c_int, c_int, c_int, array_2d_double, array_2d_double, POINTER(c_double),
                                POINTER(c_double), array_2d_double]
        else:
            sapolib=cdll.LoadLibrary("../../Sapo/libsapo_dyn_lib.so")
            sapolib.computeSapo_many.argtypes = [c_int, c_int, c_int, array_2d_double, array_2d_double, POINTER(c_double),
                                POINTER(c_double), array_2d_double,POINTER(c_float),c_int]
            #if boolean_flag:
            coeffs=poly_approx()
            deg=coeffs.shape[0]
            #    boolean_flag=False
        sapolib.computeSapo_small.restype = int
        sapolib.computeSapo_many.restype = int

        start = time.time()

        # krelu on variables in varsid
        # self.varsid = varsid
        self.k = len(cdd_hrepr[0]) - 1
        self.cdd_hrepr = cdd_hrepr
        # print("LENGTH ", len(cdd_hrepr[0]))
        # cdd_hrepr = self.get_ineqs(varsid)
        check_pt1 = time.time()

        input_cons = np.asarray(cdd_hrepr, dtype=np.double)  # this is a list of constraints
        # Ax + b >=0 with  A = input_cons[1:3][i]   b = input_cons[0][i]
        dim = input_cons.shape
        n_var = dim[1] - 1
        n_dir = dim[0] // 2
        # added call to polyfit
        if n_var == 1:
            output_cons = np.concatenate((np.tanh(input_cons[:, [0]]), input_cons[:, [1]]), axis=1)
            n_cons = dim[0]
        else:
            modelSapo = py_sapo(n_var, n_dir, input_cons)
            output_cons_temp = np.empty([dim[0], n_var + 1], dtype=np.double)
            output_cons_val = np.empty(0, dtype=np.double)

            [L, T, n_bundle] = modelSapo.LTmatrix()
            #modelSapo.offset(input_cons, dim[0])

            cL = (L.__array_interface__['data'][0]
                  + np.arange(L.shape[0]) * L.strides[0]).astype(np.uintp)
            cT = (T.__array_interface__['data'][0]
                  + np.arange(T.shape[0]) * T.strides[0]).astype(np.uintp)
            cA = (output_cons_temp.__array_interface__['data'][0]
                  + np.arange(output_cons_temp.shape[0]) * output_cons_temp.strides[0]).astype(np.uintp)


       
               
            
            if config.splitting:
                regions = modelSapo.createRegions()
                for i in range(pow(3, n_var)):
                    temp_cdd = cdd_hrepr.copy()
                    for j in range(2*n_var):
                        temp_cdd.append(regions[j + i * 2*n_var])
                    temp_cdd = cdd.Matrix(temp_cdd, number_type='fraction')
                    temp_cdd.rep_type = cdd.RepType.INEQUALITY
                    pts = cdd.Polyhedron(temp_cdd).get_generators()
                    pts_np_temp = np.array(pts, dtype=np.double)
                    
                    if len(pts_np_temp) > 0:
                        print('Region', i + 1, 'is not empty!')
                        pts_np = pts_np_temp[::, 1::]
                        if i in [0, 2, 6, 8, 18, 20, 24, 26]:
                            n_cons = 2*n_dir #pow(3, n_var) - 1
                            output_cons_val_temp = modelSapo.comput_valOutputcons(i)
                            output_cons = modelSapo.emptyoutputcons()
                            output_cons_val = np.concatenate((output_cons_val, output_cons_val_temp), axis=0)

                        else:

                            # Reshape the input constraints
                            pts_np = pts_np.transpose()
                            val = L @ pts_np
                            offp_temp = np.max(val, 1) #Lx <= b
                            offm_temp = np.max(-val, 1) #-Lx <= b

                            # Call sapo
                            coffp = offp_temp.ctypes.data_as(POINTER(c_double))
                            coffm = offm_temp.ctypes.data_as(POINTER(c_double))

    
                            if config.poly_dynamic is False:
                                n_cons = sapolib.computeSapo_small(n_var, n_dir, n_bundle, cL, cT, coffp, coffm, cA)
                            else: # add coeffs
                                c_coeffs=coeffs.ctypes.data_as(POINTER(c_float))
                                n_cons = sapolib.computeSapo_many(n_var, n_dir, n_bundle, cL, cT, coffp, coffm, cA, c_coeffs, deg)


                        # Reshape the output constraints (restrict to [-1,1]^n_var)
                        # Ax + b >= 0
                        # x_1 >= -1   x_1+1>=0
                        # x1 <= 1     -x_1+1>=0
                        # --------------------
                        # [b A] such Ax+b>=0 (se n_var=2)
                        # b0 1 0
                        # b1 0 1
                        # b2 1 1
                        # b3 1 -1
                        # b4 -1 0
                        # b5 0 -1
                        # b6 -1 -1
                        # b7 -1 1


                            if config.sanity_check:
                                if n_var == 2:
                                    output_cons_temp[[0, 1, n_dir, n_dir + 1], 0] = np.maximum(np.minimum(
                                        output_cons_temp[[0, 1, n_dir, n_dir + 1], 0], 1), -1)
                                    output_cons_temp[[2, 3, n_dir+2, n_dir+3], 0] = np.maximum(np.minimum(
                                        output_cons_temp[[2, 3, n_dir+2, n_dir+3], 0], 2), -2)
                                elif n_var == 3:
                                    output_cons_temp[[0, 1, 2, n_dir, n_dir+1, n_dir+2], 0] = np.maximum(np.minimum(
                                        output_cons_temp[[0, 1, 2, n_dir, n_dir+1, n_dir+2], 0], 1), -1)
                                    output_cons_temp[[3, 4, 5, 6, 7, 8, n_dir + 3, n_dir + 4, n_dir + 5, n_dir + 6, n_dir + 7,
                                                    n_dir + 8], 0] = np.maximum(np.minimum(
                                        output_cons_temp[[3, 4, 5, 6, 7, 8, n_dir + 3, n_dir + 4, n_dir + 5, n_dir + 6, n_dir + 7,
                                                        n_dir + 8], 0], 2),-2)
                                    output_cons_temp[
                                        [9, 10, 11, 12, n_dir + 9, n_dir + 10, n_dir + 11, n_dir + 12], 0] = np.maximum(np.minimum(
                                        output_cons_temp[[9, 10, 11, 12, n_dir + 9, n_dir + 10, n_dir + 11, n_dir + 12], 0], 3), -3)
                            else:
                                print('\nNo sanity check was performed\n')


                            # Append the bounds
                            output_cons_val = np.concatenate((output_cons_val, output_cons_temp[:, 0]), axis=0)
                            output_cons = np.copy(output_cons_temp)
           
                # Make the union of the output sets
                output_cons_val = np.reshape(output_cons_val, (-1, n_cons))
                output_cons_val = np.max(output_cons_val, 0)
                output_cons[:, 0] = output_cons_val
            else:
                # No splitting
                # Call sapo
              
                offp_temp = modelSapo.offp
                offm_temp = modelSapo.offm
                coffp = offp_temp.ctypes.data_as(POINTER(c_double))
                coffm = offm_temp.ctypes.data_as(POINTER(c_double))
                if config.poly_dynamic is False:
                    n_cons = sapolib.computeSapo_small(n_var, n_dir, n_bundle, cL, cT, coffp, coffm, cA)
                else:  # add coeffs
                    c_coeffs = coeffs.ctypes.data_as(POINTER(c_float))
                    n_cons = sapolib.computeSapo_many(n_var, n_dir, n_bundle, cL, cT, coffp, coffm, cA, c_coeffs, deg)

                #output_cons_val = np.reshape(output_cons_temp, (-1, n_cons))
                #output_cons_val = np.max(output_cons_val, 0)
                output_cons = np.copy(output_cons_temp)
               # output_cons[:, 0] = output_cons_val


        # Collect all the input-output constraints
        elaborate_input_cons = np.concatenate((input_cons, np.zeros([dim[0], n_var], dtype=np.double)), axis=1)
        elaborate_output_cons = np.concatenate((output_cons[:, [0]], np.zeros([n_cons, n_var], dtype=np.double),
                                                output_cons[:, range(n_var)]), axis=1)
        cons = np.concatenate((elaborate_input_cons, elaborate_output_cons), axis=0)

        '''
        # We get orthant points using exact precision, because it allows to guarantee soundness of the algorithm.
        cdd_hrepr = cdd.Matrix(cdd_hrepr, number_type='fraction')
        cdd_hrepr.rep_type = cdd.RepType.INEQUALITY
        pts = self.get_orthant_points(cdd_hrepr)

        df = pd.DataFrame(pts)
        df.to_csv('filename.csv', index=False)
        # Generate extremal points in the space of variables before and
        # after relu
        # HERE is the point to be changed ELE!
        pts = [([1] + row + [x if x > 0 else 0 for x in row]) for row in pts]

        adjust_constraints_to_make_sound = False
        # Floating point CDD is much faster then the precise CDD, however for some inputs it fails
        # due to numerical errors. If that is the case we fall back to using precise CDD.
        try:
            cdd_vrepr = cdd.Matrix(pts, number_type='float')
            cdd_vrepr.rep_type = cdd.RepType.GENERATOR
            # Convert back to H-repr.
            cons = cdd.Polyhedron(cdd_vrepr).get_inequalities()
            adjust_constraints_to_make_sound = True
            # I don't adjust linearities, so just setting lin_set to an empty set.
            self.lin_set = frozenset([])
        except:
            cdd_vrepr = cdd.Matrix(pts, number_type='fraction')
            cdd_vrepr.rep_type = cdd.RepType.GENERATOR
            # Convert back to H-repr.
            cons = cdd.Polyhedron(cdd_vrepr).get_inequalities()
            self.lin_set = cons.lin_set

        '''
        cons = np.asarray(cons, dtype=np.float64)

        # If floating point CDD was run, then we have to adjust constraints to make sure taht
        if 0:#adjust_constraints_to_make_sound:
            pts = np.asarray(pts, dtype=np.float64)
            cons_abs = np.abs(cons)
            pts_abs = np.abs(pts)
            cons_x_pts = np.matmul(cons, np.transpose(pts))
            cons_x_pts_err = np.matmul(cons_abs, np.transpose(pts_abs))
            # Since we use double precision number of bits to represent fraction is 52.
            # I'll use generous over-approximation by using 2^-40 as a relative error coefficient.
            rel_err = pow(2, -40)
            cons_x_pts_err *= rel_err
            cons_x_pts -= cons_x_pts_err
            for ci in range(len(cons)):
                min_val = np.min(cons_x_pts[ci, :])
                if min_val < 0:
                    cons[ci, 0] -= min_val

        # normalize constraints for numerical stability
        # more info: http://files.gurobi.com/Numerics.pdf
        #absmax = np.absolute(cons).max(axis=1)
        #self.cons = cons / absmax[:, None]

        end = time.time()

        return

    def get_orthant_points(self, cdd_hrepr):
        # Get points of polytope restricted to all possible orthtants
        pts = []
        for polarity in itertools.product([-1, 1], repeat=self.k):
            hrepr = cdd_hrepr.copy()

            # add constraints to restrict to +ve/-ve half of variables
            for i in range(self.k):
                row = [0] * (self.k + 1)
                row[1 + i] = polarity[i]
                # row corresponds to the half-space x_i>=0 if polarity[i]==+1
                hrepr.extend([row])

            # remove reduntant constraints
            hrepr.canonicalize()
            # Convert to V-repr.
            pts_new = cdd.Polyhedron(hrepr).get_generators()
            assert all(row[0] == 1 for row in pts_new)

            for row in pts_new:
                pts.append(list(row[1:]))

        return pts


def make_krelu_obj(varsid):
    return Krelu(varsid)


class Krelu_expr:
    def __init__(self, expr, varsid, bound):
        self.expr = expr
        self.varsid = varsid
        self.bound = bound


def get_ineqs_zono(varsid):
    cdd_hrepr = []

    # Get bounds on linear expressions over variables before relu
    # Order of coefficients determined by logic here
    # for coeffs in itertools.product([-2, -1, 0, 1, 2], repeat=len(varsid)):
    #     if all((abs(c) == 2 or c == 0) for c in coeffs): #ELE for redundant constraints with 2,1
    #         continue
    for coeffs in itertools.product([-2, -1, 0, 1, 2], repeat=len(varsid)):
        if len(varsid) == 1:
            if coeffs[0] == 0:
                continue
        else:
            if all((abs(c) == 2 or c == 0) for c in coeffs):  # ELE for redundant constraints with 2,1
                continue
        linexpr0 = generate_linexpr0(Krelu.offset, varsid, coeffs)
        element = elina_abstract0_assign_linexpr_array(Krelu.man, True, Krelu.element, Krelu.tdim, linexpr0, 1, None)
        bound_linexpr = elina_abstract0_bound_dimension(Krelu.man, Krelu.element, Krelu.offset + Krelu.length)
        upper_bound = bound_linexpr.contents.sup.contents.val.dbl
        cdd_hrepr.append([upper_bound] + [-c for c in coeffs])
    return cdd_hrepr


def compute_bound(constraint, lbi, ubi, varsid, j, is_lower):
    k = len(varsid)
    divisor = -constraint[j + k + 1]
    actual_bound = constraint[0] / divisor
    potential_improvement = 0
    for l in range(k):
        coeff = constraint[l + 1] / divisor
        if is_lower:
            if coeff < 0:
                actual_bound += coeff * ubi[varsid[l]]
            elif coeff > 0:
                actual_bound += coeff * lbi[varsid[l]]
        else:
            if coeff < 0:
                actual_bound += coeff * lbi[varsid[l]]
            elif coeff > 0:
                actual_bound += coeff * ubi[varsid[l]]
        potential_improvement += abs(coeff * (ubi[varsid[l]] - lbi[varsid[l]]))
        if l == j:
            continue
        coeff = constraint[l + k + 1] / divisor
        if ((is_lower and coeff < 0) or ((not is_lower) and (coeff > 0))):
            actual_bound += coeff * ubi[varsid[l]]
    return actual_bound, potential_improvement


def calculate_nnz(constraint, k):
    nnz = 0
    for i in range(k):
        if constraint[i + 1] != 0:
            nnz = nnz + 1
    return nnz


def compute_expr_bounds_from_candidates(krelu_inst, varsid, bound_expr, lbi, ubi, candidate_bounds, is_lower):
    assert not is_lower
    k = krelu_inst.k
    cons = krelu_inst.cons
    for j in range(k):
        candidate_rows = candidate_bounds[j]
        if is_lower:
            best_bound = -math.inf
        else:
            best_bound = math.inf
        best_index = -1
        for i in range(len(candidate_rows)):
            row_index = candidate_rows[i]
            actual_bound, potential_improvement = compute_bound(cons[row_index], lbi, ubi, varsid, j, is_lower)
            bound = actual_bound - potential_improvement / 2
            nnz = calculate_nnz(cons[row_index], k)
            if nnz < 2:
                continue
            if ((is_lower and bound > best_bound) or ((not is_lower) and bound < best_bound)):
                best_index = row_index
                best_bound = bound
        if best_index == -1:
            continue
        res = np.zeros(k + 1)
        best_row = cons[best_index]
        divisor = -best_row[j + k + 1]
        assert divisor > 0
        # if divisor == 0:
        #    print("ROW ",best_row)
        #    print("CONS ", cons, krelu_inst)
        #    print("CDD ", krelu_inst.cdd_hrepr)
        #    print("j ", j, "lb ", lbi[varsid[0]], lbi[varsid[1]], "ub ", ubi[varsid[0]], ubi[varsid[1]] )
        #    print("candidates ", len(candidate_rows))
        res[0] = best_row[0] / divisor
        for l in range(k):
            res[l + 1] = best_row[l + 1] / divisor
            if (l == j):
                continue
            coeff = best_row[l + k + 1] / divisor
            if ((is_lower and coeff < 0) or ((not is_lower) and (coeff > 0))):
                res[0] = res[0] + coeff * ubi[varsid[l]]
                print("res ", res, "best_row ", best_row, "j ", j)
        if varsid[j] in bound_expr.keys():
            current_bound = bound_expr[varsid[j]].bound
            if (is_lower and best_bound > current_bound) or ((not is_lower) and best_bound < current_bound):
                bound_expr[varsid[j]] = Krelu_expr(res, varsid, best_bound)
        else:
            bound_expr[varsid[j]] = Krelu_expr(res, varsid, best_bound)


def compute_expr_bounds(krelu_inst, varsid, lower_bound_expr, upper_bound_expr, lbi, ubi):
    cons = krelu_inst.cons
    nbrows = len(cons)
    k = len(varsid)
    candidate_lower_bounds = []
    candidate_upper_bounds = []
    for j in range(k):
        candidate_lower_bounds.append([])
        candidate_upper_bounds.append([])
    lin_size = len(krelu_inst.lin_set)
    new_cons = np.zeros((lin_size, 2 * k + 1), dtype=np.float64)
    lin_count = 0
    for i in range(nbrows):
        if i in krelu_inst.lin_set:
            row = cons[i]
            for j in range(2 * k + 1):
                new_cons[lin_count][j] = -row[j]
            lin_count = lin_count + 1

    krelu_inst.cons = np.vstack([cons, new_cons])
    cons = krelu_inst.cons
    nbrows = len(cons)
    for i in range(nbrows):
        row = cons[i]
        for j in range(k):
            if row[j + k + 1] < 0:
                candidate_upper_bounds[j].append(i)
            elif row[j + k + 1] > 0:
                candidate_lower_bounds[j].append(i)
    # compute_expr_bounds_from_candidates(krelu_inst, varsid, lower_bound_expr, lbi, ubi, candidate_lower_bounds, True)
    compute_expr_bounds_from_candidates(krelu_inst, varsid, upper_bound_expr, lbi, ubi, candidate_upper_bounds, False)


def get_sparse_cover_for_group_of_vars(vars):
    """Function is fast for len(vars) = 50 and becomes slow for len(vars) ~ 100."""
    K = 3
    assert len(vars) > K

    sparsed_combs = []

    for comb in itertools.combinations(vars, K):
        add = True
        for selected_comb in sparsed_combs:
            if len(set(comb).intersection(set(selected_comb))) >= K - 1:
                add = False
                break
        if add:
            sparsed_combs.append(comb)

    return sparsed_combs


def sparse_heuristic_with_cutoff(all_vars, areas):
    assert len(all_vars) == len(areas)
    K = 3
    sparse_n = config.sparse_n
    cutoff = 0.05
    print("sparse n", sparse_n)
    # Sort vars by descending area
    all_vars = sorted(all_vars, key=lambda var: -areas[var])

    vars_above_cutoff = [i for i in all_vars if areas[i] >= cutoff]

    krelu_args = []
    while len(vars_above_cutoff) > 0:
        grouplen = min(sparse_n, len(vars_above_cutoff))
        group = vars_above_cutoff[:grouplen]
        vars_above_cutoff = vars_above_cutoff[grouplen:]
        if grouplen <= K:
            krelu_args.append(group)
        else:
            group_args = get_sparse_cover_for_group_of_vars(group)

            for arg in group_args:
                krelu_args.append(arg)

    # Also just apply 1-relu for every var.
    for var in all_vars:
        krelu_args.append([var])

    return krelu_args


def encode_kactivation_cons(nn, man, element, offset, layerno, length, lbi, ubi, relu_groups, need_pop, domain,
                            activation_type):
    import deepzono_nodes as dn
    if (need_pop):
        relu_groups.pop()

    last_conv = -1
    is_conv = False
    for i in range(nn.numlayer):
        if nn.layertypes[i] == 'Conv':
            last_conv = i
            is_conv = True

    lbi = np.asarray(lbi, dtype=np.double)
    ubi = np.asarray(ubi, dtype=np.double)

    candidate_vars = [i for i in range(length) if lbi[i] < 0 and ubi[i] > 0]
    candidate_vars_areas = {var: -lbi[var] * ubi[var] for var in candidate_vars}
    # Sort vars by descending area
    candidate_vars = sorted(candidate_vars, key=lambda var: -candidate_vars_areas[var])

    # Use sparse heuristic to select args (uncomment to use)
    krelu_args = sparse_heuristic_with_cutoff(candidate_vars, candidate_vars_areas)

    relucons = []
    # print("UBI ",ubi)
    tdim = ElinaDim(offset + length)
    if domain == 'refinezono':
        element = dn.add_dimensions(man, element, offset + length, 1)

    # krelu_args = []
    # if config.dyn_krelu and candidate_vars:
    #    limit3relucalls = 500
    #    firstk = math.sqrt(6*limit3relucalls/len(candidate_vars))
    #    firstk = int(min(firstk, len(candidate_vars)))
    #    if is_conv and layerno < last_conv:
    #        firstk = 1
    #    else:
    #        firstk = 5#int(max(1,firstk))
    #    print("firstk ",firstk)
    #    if firstk>3:
    #        while candidate_vars:
    #            headlen = min(firstk, len(candidate_vars))
    #            head = candidate_vars[:headlen]
    #            candidate_vars = candidate_vars[headlen:]
    #            if len(head)<=3:
    #               krelu_args.append(head)
    #            else:
    #                for arg in itertools.combinations(head, 3):
    #                    krelu_args.append(arg)

    # klist = ([3] if (config.use_3relu) else []) + ([2] if (config.use_2relu) else []) + [1]
    # for k in klist:
    #    while len(candidate_vars) >= k:
    #        krelu_args.append(candidate_vars[:k])
    #        candidate_vars = candidate_vars[k:]
    Krelu.man = man
    Krelu.element = element
    Krelu.tdim = tdim
    Krelu.length = length
    Krelu.layerno = layerno
    Krelu.offset = offset
    Krelu.domain = domain

    start = time.time()
    if domain == 'refinezono':
        with multiprocessing.Pool(config.numproc) as pool:
            cdd_hrepr_array = pool.map(get_ineqs_zono, krelu_args)
    else:
        #    krelu_results = []
        total_size = 0
        listCoeff = [-2, -1, 0, 1, 2]
        for varsid in krelu_args:
            if len(varsid) == 1:
                size = len(listCoeff) ** len(varsid) - 1
            else:
                size = len(listCoeff) ** len(varsid) - (3 ** len(varsid))
            total_size = total_size + size

        linexpr0 = elina_linexpr0_array_alloc(total_size)
        # HERE we should define our expressions
        i = -1
        for varsid in krelu_args:
            for coeffs in itertools.product(listCoeff, repeat=len(varsid)):
                i = i + 1
                if len(varsid) == 1:
                    if coeffs[0] == 0:
                        i = i - 1
                        continue
                else:
                    if all((abs(c) == 2 or c == 0) for c in coeffs): #ELE for redundant constraints with 2,1
                        i = i - 1
                        continue

                linexpr0[i] = generate_linexpr0(offset, varsid, coeffs)
        upper_bound = get_upper_bound_for_linexpr0(man, element, linexpr0, total_size, layerno)
        i = -1
        cdd_hrepr_array = []
        bound_val = []
        for varsid in krelu_args:
            cdd_hrepr = []
            for coeffs in itertools.product(listCoeff, repeat=len(varsid)):
                i = i + 1
                if len(varsid) == 1:
                    if coeffs[0] == 0:
                        i = i - 1
                        continue
                else:
                    if all((abs(c) == 2 or c == 0) for c in coeffs):  # ELE for redundant constraints with 2,1
                        i = i - 1
                        continue
                cdd_hrepr.append([upper_bound[i]] + [-c for c in coeffs])
                bound_val.append(upper_bound[i])  # ELE
            #     print("UPPER BOUND ", upper_bound[i], "COEFF ", coeffs)
            # if len(varsid)>1:
            #     print("LB ", lbi[varsid[0]],lbi[varsid[1]], "UB ", ubi[varsid[0]], ubi[varsid[1]])
            cdd_hrepr_array.append(cdd_hrepr)

    with multiprocessing.Pool(config.numproc) as pool:  # here is entering in 'get_orthant_points'
        krelu_results = pool.map(make_krelu_obj, cdd_hrepr_array)#, lbi, ubi)

    #        krelu_results.append(make_krelu_obj(krelu_args[i]))
    # bound_expr_list = []
    gid = 0
    lower_bound_expr = {}
    upper_bound_expr = {}
    for krelu_inst in krelu_results:
        varsid = krelu_args[gid]
        krelu_inst.varsid = varsid
        # For now disabling since in the experiments updating expression bounds makes results worse.
        # compute_expr_bounds(krelu_inst, varsid, lower_bound_expr, upper_bound_expr, lbi, ubi)
        # print("VARSID ",varsid)
        # bound_expr_list.append(Krelu_expr(lower_bound_expr, upper_bound_expr, varsid))
        relucons.append(krelu_inst)
        gid = gid + 1
    end = time.time()

    if config.debug:
        print('krelu time spent: ' + str(end - start))
    if domain == 'refinezono':
        element = dn.remove_dimensions(man, element, offset + length, 1)

    relu_groups.append(relucons)

    return lower_bound_expr, upper_bound_expr
