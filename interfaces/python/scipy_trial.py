import scipy.sparse as sc
import numpy as np

row = np.array([0, 2, 2, 0, 1, 2])
col = np.array([0, 0, 1, 2, 2, 2])
data = np.array([1, 2, 3, 4, 5, 6])
Q = sc.csc_matrix((data, (row, col)), shape=(3, 3))
Q_sym = (Q+Q.transpose())/2
Q_array = Q.toarray()

print(Q_sym.indptr)
print(Q_sym.toarray())

(n,m) = Q_sym.shape
print(n)

q = np.array([0, 2, 3])
print(len(q))