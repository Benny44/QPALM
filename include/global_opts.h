/**
 * @file global_opts.h
 * @author Ben Hermans
 * @brief Custom memory allocation, print and utility functions, and data types for floats and ints.
 * @details Memory allocation and print functions depend on whether the code is compiled as a standalone
 * library or with matlab or python. The data types used for floating point numbers and integer numbers
 * can be changed here as well. Finally, some customized operations (max, min, mod and abs) are included
 * as well.
 */

#ifndef GLOBAL_OPTS_H
# define GLOBAL_OPTS_H

# ifdef __cplusplus
extern "C" {
# endif /* ifdef __cplusplus */

#include "cholmod.h"

/**
 * @name CHOLMOD data types and macros 
 * @{
 */
#ifdef DLONG
#define Real double
#define Int SuiteSparse_long
#define Int_max SuiteSparse_long_max
#define CHOLMOD(name) cholmod_l_ ## name
//#define LONG
//#define DOUBLE
#define ITYPE CHOLMOD_LONG
#define DTYPE CHOLMOD_DOUBLE
#define ID SuiteSparse_long_id

#else

#ifndef DINT
#define DINT
#endif
//#define INT
//#define DOUBLE

#define Real double
#define Int int
#define Int_max INT_MAX
#define CHOLMOD(name) cholmod_ ## name
#define ITYPE CHOLMOD_INT
#define DTYPE CHOLMOD_DOUBLE
#define ID "%d"

/* GPU acceleration is not available for the int version of CHOLMOD */
#undef GPU_BLAS

#endif

/**
 * @}
 */

/** 
 * @name Data customizations
 * @{
 */
#  include <stdlib.h>

typedef Real c_float; /**< doubles for numerical values  */
typedef Int c_int; /**< long or int for indices */

/**
 * @}
 */

/**
 * @name Custom memory allocation (e.g. matlab/python)
 * @{
 */
#  ifdef MATLAB
    #   include "mex.h"
static void* c_calloc(size_t num, size_t size) {
  void *m = mxCalloc(num, size);
  mexMakeMemoryPersistent(m);
  return m;
}

static void* c_malloc(size_t size) {
  void *m = mxMalloc(size);

  mexMakeMemoryPersistent(m);
  return m;
}

static void* c_realloc(void *ptr, size_t size) {
  void *m = mxRealloc(ptr, size);

  mexMakeMemoryPersistent(m);
  return m;
}

    #   define c_free mxFree
#  elif defined PYTHON

// Define memory allocation for python. Note that in Python 2 memory manager
// Calloc is not implemented
    #   include <Python.h>
    #   define c_malloc PyMem_Malloc
    #   if PY_MAJOR_VERSION >= 3
    #    define c_calloc PyMem_Calloc
    #   else  /* if PY_MAJOR_VERSION >= 3 */
static void* c_calloc(size_t num, size_t size) {
  void *m = PyMem_Malloc(num * size);

  memset(m, 0, num * size);
  return m;
}

    #   endif /* if PY_MAJOR_VERSION >= 3 */

// #define c_calloc(n,s) ({
//         void * p_calloc = c_malloc((n)*(s));
//         memset(p_calloc, 0, (n)*(s));
//         p_calloc;
//     })
    #   define c_free PyMem_Free
    #   define c_realloc PyMem_Realloc
#  else  /* if not MATLAB of Python */
    #   define c_malloc malloc    /**< custom malloc */
    #   define c_calloc calloc    /**< custom calloc */
    #   define c_free free        /**< custom free */
    #   define c_realloc realloc  /**< custom realloc */

#  endif /* ifdef MATLAB */

/**
 * @}
 */


/* PRINTING */
# ifdef PRINTING
#  include <stdio.h>
#  include <string.h>

#  ifdef MATLAB
#   define c_print mexPrintf

#  elif defined PYTHON
#   include <Python.h>
#   define c_print PySys_WriteStdout
#  elif defined R_LANG
#   include <R_ext/Print.h>
#   define c_print Rprintf
#  else  /* ifdef MATLAB */
#   define c_print printf
#  endif /* ifdef MATLAB */

// Print error macro
#  define c_eprint(...) c_print("ERROR in %s: ", __FUNCTION__); c_print( \
    __VA_ARGS__); c_print("\n");

# endif /* ifdef PRINTING */


/**
 * @name Custom operations
 * @{
 */
# ifndef c_absval
#  define c_absval(x) (((x) < 0) ? -(x) : (x)) /**< absolute value */
# endif /* ifndef c_absval */

# ifndef c_max
#  define c_max(a, b) (((a) > (b)) ? (a) : (b)) /**< maximum of two values */
# endif /* ifndef c_max */

# ifndef c_min
#  define c_min(a, b) (((a) < (b)) ? (a) : (b)) /**< minimum of two values */
# endif /* ifndef c_min */

# ifndef mod
#  define mod(a,b) ((((a)%(b))+(b))%(b)) /**< modulo operation (positive result for all values) */
#endif

#include <math.h>
#  define c_sqrt sqrt /**< square root */

/** @} */

# ifdef __cplusplus
}
# endif /* ifdef __cplusplus */

#endif /* ifndef GLOBAL_OPTS_H */