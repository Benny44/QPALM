from ctypes import *
import platform
import numpy as np
import scipy as sc


class QPALMSettings(Structure):
    _fields_ = [("max_iter", c_int),
                ("eps_abs", c_double),
                ("eps_rel", c_double),
                ("eps_abs_in", c_double),
                ("eps_rel_in", c_double),
                ("rho", c_double),
                ("eps_prim_inf", c_double),
                ("eps_dual_inf", c_double),
                ("theta", c_double),
                ("delta", c_double),
                ("tau_init", c_double),
                ("proximal", c_int),
                ("gamma_init", c_double),
                ("gamma_upd", c_double),
                ("gamma_max", c_double),
                ("scaling", c_int),
                ("nonconvex", c_int),
                ("verbose", c_int),
                ("warm_start", c_int),
                ("time_limit", c_double)
                ]

QPALMSettings_pointer = POINTER(QPALMSettings)

class cholmod_sparse(Structure):
    _fields_ = [("nrow", c_uint),
                ("ncol", c_uint),
                ("nzmax", c_uint),
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

class cholmod_common(Structure):
    _fields_ = [("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ("", ),
                ]


class QPALMData(Structure):
    _fields_ = [("n", c_uint),
                ("m", c_uint),
                ("Q", cholmod_sparse_pointer),
                ("A", cholmod_sparse_pointer),
                ("q", POINTER(c_double)),
                ("c", c_double),
                ("bmin", POINTER(c_double)),
                ("bmax", POINTER(c_double))
                ]
QPALMData_pointer = POINTER(QPALMData)



class Qpalm:
    """
    Wrapper class for the python interface to QPALM
    """
    def __init__(self):
        """
        Construct the wrapper class and load the dynamic library.
        """
        self._load_library()
        self._set_restypes()
        self._settings = self.python_interface.qpalm_malloc_settings()

    #def __del__(self):
        #self.python_interface.qpalm_cleanup(work)
    def set_default_settings(self):
        self.python_interface.qpalm_set_default_settings(self._settings)
        
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

        c_double_p = POINTER(c_double)
        c_int_p = POINTER(c_int)

        self._data[0].n = n
        self._data[0].m = m
        self._data[0].q = q.ctypes.as_data(c_double_p)
        self._data[0].bmin = bmin.ctypes.as_data(c_double_p)
        self._data[0].bmax = bmax.ctypes.as_data(c_double_p)

        self._data[0].A[0].nrow = m
        self._data[0].A[0].ncol = n
        Ap = A.indptr
        Ai = A.indices
        self._data[0].A[0].p = Ap.ctypes.as_data(c_int_p)
        self._data[0].A[0].i = Ai.ctypes.as_data(c_int_p)
        self._data[0].A[0].nzmax = Ap[n]
        self._data[0].A[0].packed = 1
        self._data[0].A[0].sorted = 1
        self._data[0].A[0].nz = 0 #NULL
        self._data[0].A[0].itype = 2 #CHOLMOD_LONG 
        self._data[0].A[0].dtype = 0 #CHOLMOD_DOUBLE
        self._data[0].A[0].stype = 0 #Unsymmetric
        Ax = A.data
        self._data[0].A[0].x = Ax.ctypes.as_data(c_double_p)
        self._data[0].A[0].xtype = 1 #CHOLMOD_REAL

        self._data[0].Q[0].nrow = n
        self._data[0].Q[0].ncol = n
        Qp = Q.indptr
        Qi = Q.indices
        self._data[0].Q[0].p = Qp.ctypes.as_data(c_int_p)
        self._data[0].Q[0].i = Qi.ctypes.as_data(c_int_p)
        self._data[0].Q[0].nzmax = Qp[n]
        self._data[0].Q[0].packed = 1
        self._data[0].Q[0].sorted = 1
        self._data[0].Q[0].nz = 0 #NULL
        self._data[0].Q[0].itype = 2 #CHOLMOD_LONG 
        self._data[0].Q[0].dtype = 0 #CHOLMOD_DOUBLE
        self._data[0].Q[0].stype = -1 #Lower symmetric
        Qx = Q.data
        self._data[0].Q[0].x = Qx.ctypes.as_data(c_double_p)
        self._data[0].Q[0].xtype = 1 #CHOLMOD_REAL

    def _allocate_work(self):

        work = self.python_interface.qpalm_setup()


    def _load_library(self):
        """
        Load the dynamic QPALM library.
        """
        try:
            if (platform.system() == 'Linux'):
                print("OS is Linux")      
                self.python_interface = CDLL("../build/lib/" + "libqpalm.so")
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

if __name__== '__main__':
    qpalm = Qpalm()
    qpalm.set_default_settings()
    qpalm._settings[0].max_iter = 4
    print(qpalm._settings[0].max_iter)