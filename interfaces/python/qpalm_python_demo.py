import qpalm as qp
import numpy as np
import scipy as sc
import scipy.sparse as sp


solver = qp.Qpalm()

row = np.array([0, 0, 1, 1])
col = np.array([0, 1, 0, 1])
data = np.array([1, -1, -1, 2])
Q = sp.csc_matrix((data, (row, col)), shape=(3, 3))

q = np.array([-2, -6, 1])
bmin = np.array([0.5, -10, -10, -10])
bmax = np.array([0.5, 10, 10, 10])

row = np.array([0, 1, 0, 2, 0, 3])
col = np.array([0, 0, 1, 1, 2, 2])
data = np.array([1, 1, 1, 1, 1, 1])
A = sp.csc_matrix((data, (row, col)), shape=(4, 3))

solver.set_data(Q=Q, A=A, q=q, bmin=bmin, bmax=bmax)
solver._allocate_work()
solver._solve()
sol_x = solver._work.contents.solution.contents.x
tol = 1e-5
assert(abs(sol_x[0] - 5.5) < tol)
assert(abs(sol_x[1] - 5.0) < tol)
assert(abs(sol_x[2] - (-10)) < tol)