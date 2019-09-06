/**
 * @file constants.h
 * @author Ben Hermans
 * @brief Macros used in QPALM
 * @details This file contains the macros that are used as default settings and to set the solver status.
 * @warning Verbose is not implemented as of yet.
 */
#ifndef CONSTANTS_H
#define CONSTANTS_H

#ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

/**
 * @name Booleans
 * @{
 */
#define TRUE 1
#define FALSE 0
/**
 * @}
 */


/** 
 * @name Solver status
 * @{
 */

# define QPALM_SOLVED (1)             /**< status to indicate the problem is solved to optimality given the specified tolerances */
# define QPALM_MAX_ITER_REACHED (-2)  /**< status to indicate termination due to reaching the maximum number of iterations  */
# define QPALM_PRIMAL_INFEASIBLE (-3) /**< status to indicate the problem is primal infeasible  */
# define QPALM_DUAL_INFEASIBLE (-4)   /**< status to indicate the problem is dual infeasible  */
# define QPALM_UNSOLVED (-10)         /**< status to indicate the problem is unsolved. Only setup function has been called */
# define QPALM_ERROR (0)              /**< status to indicate an error has occured (this error should automatically be printed) */

/**
 * @}
 */

/** 
 * @name Solver parameters and settings
 * @{
 */

/**********************************
* Solver Parameters and Settings *
**********************************/

# define MAX_ITER (10000)       /**< default maximum number of iterations */
# define INNER_MAX_ITER (100)   /**< default maximum number of iterations per subproblem */
# define EPS_ABS (1e-4)         /**< default absolute convergence tolerance */
# define EPS_REL (1e-4)         /**< default relative convergence tolerance */
# define EPS_ABS_IN (1)         /**< default intermediate absolute convergence tolerance */
# define EPS_REL_IN (1)         /**< default intermediate relative convergence tolerance */
# define RHO (0.1)              /**< default tolerance scaling factor */
# define EPS_PRIM_INF (1e-5)    /**< default primal infeasibility tolerance */
# define EPS_DUAL_INF (1e-5)    /**< default dual infeasibility tolerance */
# define THETA (0.25)           /**< default penalty update criterion parameter */
# define DELTA (100)             /**< default penalty update factor */
# define TAU_INIT (1)           /**< default initial stepsize in backtracking */

# define PROXIMAL (TRUE)        /**< default use of proximal method of multipliers */
# define GAMMA_INIT (1E1)       /**< default initial proximal penalty parameter */
# define GAMMA_UPD (10)         /**< default proximal penalty update factor */
# define GAMMA_MAX (1E7)        /**< default proximal penalty cap */

# define SCALING (2)            /**< default number of scaling iterations */
# define MIN_SCALING (1e-12)    /**< minimum scaling value *////< Minimum scaling value
# define MAX_SCALING (1e+04)    /**< maximum scaling value *////< Maximum scaling value

# define NONCONVEX (FALSE)       /**< default use of nonconvex adjustments */
# define WARM_START (FALSE)     /**< default warm start setting */
# define VERBOSE (TRUE)         /**< default write out progress setting */

# define MAX_RANK_UPDATE 160 /**< maximum rank for the sparse factorization update */

# ifndef QPALM_NULL
#  define QPALM_NULL 0 /**< NULL, if something goes wrong during setup, the workspace pointer is set to this */
# endif /* ifndef QPALM_NULL */

# ifndef QPALM_NAN
#  define QPALM_NAN ((c_float)0x7fc00000UL)  /**< not a number, used for the solution if the problem is primal or dual infeasible */
# endif /* ifndef QPALM_NAN */

# ifndef QPALM_INFTY
#  define QPALM_INFTY ((c_float)1e20)  /**< infinity, used to indicate one-sided constraints */
# endif /* ifndef QPALM_INFTY */

/**
 * @}
 */


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CONSTANTS_H