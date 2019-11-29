from ctypes import *
import platform
import numpy as np
import scipy as sc
import scipy.sparse as sp
import os

class QPALMSettings(Structure):
    _fields_ = [("max_iter", c_long),
                ("inner_max_iter", c_long),
                ("eps_abs", c_double),
                ("eps_rel", c_double),
                ("eps_abs_in", c_double),
                ("eps_rel_in", c_double),
                ("rho", c_double),
                ("eps_prim_inf", c_double),
                ("eps_dual_inf", c_double),
                ("theta", c_double),
                ("delta", c_double),
                ("sigma_max", c_double),
                ("proximal", c_long),
                ("gamma_init", c_double),
                ("gamma_upd", c_double),
                ("gamma_max", c_double),
                ("scaling", c_long),
                ("nonconvex", c_long),
                ("verbose", c_long),
                ("print_iter", c_long),
                ("warm_start", c_long),
                ("reset_newton_iter", c_long),
                ("enable_dual_termination", c_long),
                ("dual_objective_limit", c_double),
                ("time_limit", c_double)
                ]

QPALMSettings_pointer = POINTER(QPALMSettings)

c_long_p = POINTER(c_long)
c_double_p = POINTER(c_double)

class cholmod_sparse(Structure):
    _fields_ = [("nrow", c_ulong),
                ("ncol", c_ulong),
                ("nzmax", c_ulong),
                ("p", c_void_p),
                ("i", c_void_p),
                ("nz", c_void_p),
                ("x", c_void_p),
                ("z", c_void_p),
                ("stype", c_int),
                ("itype", c_int),
                ("xtype", c_int),
                ("dtype", c_int),
                ("sorted", c_int),
                ("packed", c_int)
                ]

cholmod_sparse_pointer = POINTER(cholmod_sparse)

class QPALMData(Structure):
    _fields_ = [("n", c_ulong),
                ("m", c_ulong),
                ("Q", cholmod_sparse_pointer),
                ("A", cholmod_sparse_pointer),
                ("q", POINTER(c_double)),
                ("c", c_double),
                ("bmin", POINTER(c_double)),
                ("bmax", POINTER(c_double))
                ]
QPALMData_pointer = POINTER(QPALMData)


class QPALMInfo(Structure):
    _fields_ = [("iter", c_long),
                ("iter_out", c_long),
                ("status", c_char*32),
                ("status_val", c_long),
                ("pri_res_norm", c_double),
                ("dua_res_norm", c_double),
                ("dua2_res_norm", c_double),
                ("objective", c_double),
                ("dual_objective", c_double),
                ("setup_time", c_double),
                ("solve_time", c_double),
                ("run_time", c_double)
                ]

QPALMInfo_pointer = POINTER(QPALMInfo)

class array_element(Structure):
    _fields_ = [("x", c_double),
                ("i", c_long)]

class QPALMSolution(Structure):
    _fields_ = [("x", POINTER(c_double)),
                ("y", POINTER(c_double))]


class QPALMWork(Structure):
    _fields_ = [("data", QPALMData_pointer),
                ("x", POINTER(c_double)),
                ("y", POINTER(c_double)),
                ("Ax", POINTER(c_double)),
                ("Qx", POINTER(c_double)),
                ("Aty", POINTER(c_double)),
                ("x_prev", POINTER(c_double)),
                ("initialized", c_long),
                ("temp_m", POINTER(c_double)),
                ("temp_n", POINTER(c_double)),
                ("sigma", POINTER(c_double)),
                ("sqrt_sigma_max", c_double),
                ("nb_sigma_changed", c_long),
                ("gamma", c_double),
                ("gamma_maxed", c_long),
                ("Axys", POINTER(c_double)),
                ("z", POINTER(c_double)),
                ("pri_res", POINTER(c_double)),
                ("pri_res_in", POINTER(c_double)),
                ("yh", POINTER(c_double)),
                ("Atyh", POINTER(c_double)),
                ("df", POINTER(c_double)),
                ("x0", POINTER(c_double)),
                ("xx0", POINTER(c_double)),
                ("dphi", POINTER(c_double)),
                ("neg_dphi", POINTER(c_double)),
                ("dphi_prev", POINTER(c_double)),
                ("d", POINTER(c_double)),
                ("tau", c_double),
                ("Qd", POINTER(c_double)),
                ("Ad", POINTER(c_double)),
                ("sqrt_sigma", POINTER(c_double)),
                ("sqrt_delta", c_double),
                ("eta", c_double),
                ("beta", c_double),
                ("delta", POINTER(c_double)),
                ("alpha", POINTER(c_double)),
                ("temp_2m", POINTER(c_double)),
                ("delta2", POINTER(c_double)),
                ("delta_alpha", POINTER(c_double)),
                ("s", POINTER(array_element)),
                ("index_L", POINTER(c_long)),
                ("index_P", POINTER(c_long)),
                ("index_J", POINTER(c_long)),
                ("eps_pri", c_double),
                ("eps_dua", c_double),
                ("eps_dua_in", c_double),
                ("eps_abs_in", c_double),
                ("eps_rel_in", c_double),
                ("delta_y", POINTER(c_double)),
                ("Atdelta_y", POINTER(c_double)),
                ("delta_x", POINTER(c_double)),
                ("Qdelta_x", POINTER(c_double)),
                ("Adelta_x", POINTER(c_double)),
                ("D_temp", POINTER(c_double)),
                ("E_temp", POINTER(c_double)),
                ("chol", c_void_p),
                ("settings", QPALMSettings_pointer),
                ("scaling", c_void_p),
                ("solution", POINTER(QPALMSolution)),
                ("info", QPALMInfo_pointer),
                ("timer", c_void_p)
                ]

QPALMWork_pointer = POINTER(QPALMWork)


class Qpalm:
    """
    Wrapper class for the python interface to QPALM
    """
    def __init__(self):
        """
        Construct the wrapper class and load the dynamic library.
        """
        self._work = None
        self._load_library()
        self._set_restypes()
        self._settings = self.python_interface.qpalm_malloc_settings()
        self.python_interface.qpalm_set_default_settings(self._settings)
    
    def __del__(self): #TODO free the data and settings
        self.python_interface.qpalm_cleanup(self._work)
        self.python_interface.qpalm_free_settings(self._settings)
        self.python_interface.qpalm_free_data(self._data)

    # def set_default_settings(self):
    #     self.python_interface.qpalm_set_default_settings(self._settings)
        
    def set_data(self, Q, A, q, bmin, bmax):
        """
        Convert the data to QPALMData structure.
        Parameters
        ---------
        Q : Quadratic part of the cost (scipy.csc_matrix)
        A : Constraint matrix (scipy.csc_matrix)
        q : Linear part of the cost (numpy.array)
        bmin : Lower bounds of the constraints (numpy.array)
        bmax : Upper bounds of the constraints (numpy.array)
        """
        self._data = self.python_interface.qpalm_malloc_data()
        
        self._data.contents.c = 0

        (n,m) = Q.shape
        if n != m :
            print("ERROR: Q is not a square matrix")
        if len(q) != n :
            print("ERROR: q is not the right length")  
        (m,nA) = A.shape 
        if m != 0 and n != nA :
            print("ERROR: A is not the right size")
        if len(bmin) != m :
            print("ERROR: bmin is not the right length")
        if len(bmax) != m :
            print("ERROR: bmax is not the right length")            

        self._data[0].n = n
        self._data[0].m = m
        q = q.astype(np.float64)
        self._data[0].q = q.ctypes.data_as(c_double_p)
        bmin = bmin.astype(np.float64)
        self._data[0].bmin = bmin.ctypes.data_as(c_double_p)
        bmax = bmax.astype(np.float64)
        self._data[0].bmax = bmax.ctypes.data_as(c_double_p)

        #Make Q symmetric
        Q = (Q+Q.transpose())/2

        self._data[0].A = self.python_interface.python_allocate_cholmod_sparse(m, n, A.nnz)
        self._data[0].Q = self.python_interface.python_allocate_cholmod_sparse(n, n, Q.nnz)

        Ap = A.indptr
        Ap = Ap.astype(np.int64)
        Ai = A.indices
        Ai = Ai.astype(np.int64)

        self._data[0].A[0].p = Ap.ctypes.data_as(c_void_p)
        self._data[0].A[0].i = Ai.ctypes.data_as(c_void_p)

        self._data[0].A[0].stype = 0 #Unsymmetric
        Ax = A.data
        Ax = Ax.astype(np.float64)
        self._data[0].A[0].x = Ax.ctypes.data_as(c_void_p)

        Qp = Q.indptr
        Qp = Qp.astype(np.int64)
        Qi = Q.indices
        Qi = Qi.astype(np.int64)
        self._data[0].Q[0].p = Qp.ctypes.data_as(c_void_p)
        self._data[0].Q[0].i = Qi.ctypes.data_as(c_void_p)

        self._data[0].Q[0].stype = -1 #Lower symmetric
        Qx = Q.data
        Qx = Qx.astype(np.float64)
        self._data[0].Q[0].x = Qx.ctypes.data_as(c_void_p)


    def _allocate_work(self):

        self._work = self.python_interface.qpalm_setup(self._data, self._settings)

    def _solve(self):

        self.python_interface.qpalm_solve(self._work)

    def _load_library(self):
        """
        Load the dynamic QPALM library.
        """
        lib_dir = os.path.dirname(os.path.abspath(__file__)) 
        try:
            if (platform.system() == 'Linux'):
                print("OS is Linux")      
                lib_dir = lib_dir + "/../../build/lib/"
                self.python_interface = CDLL(lib_dir + "libqpalm.so")
            elif (platform.system() == 'Windows'):
                print("OS is Windows")
            elif (platform.system() == 'Darwin'):
                print("OS is MacOS")
            else:
                print("ERROR: could not detect OS, using Linux")
        except:
            print("Failed to load dynamic library")

    def _set_restypes(self):
        """
        Set the return types for the relavent interface functions.
        """
        self.python_interface.qpalm_malloc_settings.restype = QPALMSettings_pointer
        self.python_interface.qpalm_malloc_data.restype = QPALMData_pointer
        self.python_interface.qpalm_setup.restype = QPALMWork_pointer
        self.python_interface.qpalm_solve.restype = None
        self.python_interface.qpalm_cleanup.restype = None
        self.python_interface.python_allocate_cholmod_sparse.restype = cholmod_sparse_pointer
        self.python_interface.qpalm_free_settings.restype = None
        self.python_interface.qpalm_free_data.restype = None


if __name__== '__main__':
    qpalm = Qpalm()

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

    qpalm.set_data(Q=Q, A=A, q=q, bmin=bmin, bmax=bmax)
    qpalm._allocate_work()
    qpalm._solve()
    sol_x = qpalm._work.contents.solution.contents.x
    tol = 1e-5
    assert(abs(sol_x[0] - 5.5) < tol)
    assert(abs(sol_x[1] - 5.0) < tol)
    assert(abs(sol_x[2] - (-10)) < tol)
    