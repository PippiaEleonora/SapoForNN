import numpy as np
from ctypes import cdll
from ctypes import *
import numpy.ctypeslib

array_2d_double = np.ctypeslib.ndpointer(dtype=np.uintp, ndim=1, flags='C') 

sapolib = cdll.LoadLibrary("./libsapo_dyn_lib.so")
sapolib.computeSapo.argtypes = [c_int, c_int, c_int, array_2d_double, array_2d_double, POINTER(c_double), POINTER(c_double), array_2d_double,POINTER(c_float),c_int]
sapolib.computeSapo.restype = int 

#do not work in 1D
#n_var = 1
#n_dir = 1
#n_bundle = 1
#
#T = np.array([[0]], dtype=np.double)
#offp = np.array([1], dtype=np.double)
#offm = np.array([1], dtype=np.double)

n_var = 3
n_dir = 6
n_bundle = 2

A = np.empty([100, n_var+1], dtype=np.double)
L = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1]], dtype=np.double)
T = np.array([[0, 1, 2], [3, 4, 5]], dtype=np.double)
offp = np.array([1, 1, 3, 1, 2, 3], dtype=np.double)
offm = np.array([1, 0, -0.5, 0.5, -0.5, -1], dtype=np.double)
coeffs=np.array([0,1,5,0,3],dtype=np.double)

coffp = offp.ctypes.data_as(POINTER(c_double))
coffm = offm.ctypes.data_as(POINTER(c_double))
c_coeffs = coeffs.ctypes.data_as(POINTER(c_float))
deg=len(coeffs)
cL = (L.__array_interface__['data'][0] 
      + np.arange(L.shape[0])*L.strides[0]).astype(np.uintp)
cT = (T.__array_interface__['data'][0] 
      + np.arange(T.shape[0])*T.strides[0]).astype(np.uintp)
cA = (A.__array_interface__['data'][0] 
      + np.arange(A.shape[0])*A.strides[0]).astype(np.uintp)

n_cons = sapolib.computeSapo(n_var, n_dir, n_bundle, cL, cT, coffp, coffm, cA,c_coeffs,deg)

print(A[0:n_cons][:])
