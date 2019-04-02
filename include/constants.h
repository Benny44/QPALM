#ifndef CONSTANTS_H
# define CONSTANTS_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

/******************
* Solver Status  *
******************/

# define QPALM_SOLVED (1)
# define QPALM_MAX_ITER_REACHED (-2)
# define QPALM_PRIMAL_INFEASIBLE (-3) /* primal infeasible  */
# define QPALM_DUAL_INFEASIBLE (-4)   /* dual infeasible */
# define QPALM_NON_CVX (-7)           /* problem non convex */
# define QPALM_UNSOLVED (-10) /* Unsolved. Only setup function has been called */


/**********************************
* Solver Parameters and Settings *
**********************************/

# define MAX_ITER (10000)
# define EPS_ABS (1E-4)
# define EPS_REL (1E-4)
# define EPS_ABS_IN (1)
# define EPS_REL_IN (1)
# define RHO (0.1)
# define EPS_PRIM_INF (1E-4)
# define EPS_DUAL_INF (1E-4)
# define THETA (0.25)
# define DELTA (2)
# define TAU_INIT (1)

# define MEMORY (10)

# define PROXIMAL (1)
# define GAMMA (1E4)
# define GAMMA_UPD (10)
# define GAMMA_MAX (1E8)

# define SCALING (10)
# define MIN_SCALING (1e-08) ///< Minimum scaling value
# define MAX_SCALING (1e+04) ///< Maximum scaling value

# define WARM_START (0)
# define VERBOSE (1)


# ifndef QPALM_NULL
#  define QPALM_NULL 0
# endif /* ifndef QPALM_NULL */

# ifndef QPALM_NAN
#  define QPALM_NAN ((c_float)0x7fc00000UL) // Not a Number
# endif /* ifndef QPALM_NAN */

# ifndef QPALM_INFTY
#  define QPALM_INFTY ((c_float)1e20) // Infinity
# endif /* ifndef QPALM_INFTY */




# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CONSTANTS_H