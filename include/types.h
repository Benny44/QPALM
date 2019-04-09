#ifndef QPALM_TYPES_H
# define QPALM_TYPES_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "global_opts.h"
#include "cholmod.h"

/**
 * Array to sort in linesearch
 */
typedef struct array_element  {
  c_float x; ///< Array element
  c_int   i; ///< Index
} array_element;


/******************
* Internal types *
******************/

/**
 *  Matrix in compressed-column or triplet form
 */
typedef struct {
  c_int    nzmax; ///< maximum number of entries.
  c_int    m;     ///< number of rows
  c_int    n;     ///< number of columns
  c_int   *p;     ///< column pointers (size n+1) (col indices (size nzmax)
                  // start from 0 when using triplet format (direct KKT matrix
                  // formation))
  c_int   *i;     ///< row indices, size nzmax starting from 0
  c_float *x;     ///< numerical values, size nzmax
  c_int    nz;    ///< # of entries in triplet matrix, -1 for csc
} csc;


/**
 * Solution structure
 */
typedef struct {
  c_float *x; ///< Primal solution
  c_float *y; ///< Dual solution$
} QPALMSolution;

/**
 * QPALM Timer for statistics
 */
typedef struct QPALM_TIMER QPALMTimer;

/**
 * Problem scaling matrices stored as vectors
 */
typedef struct {
  c_float *D;    ///< primal variable scaling
  c_float *Dinv; ///< primal variable rescaling
  c_float *E;    ///< dual variable scaling
  c_float *Einv; ///< dual variable rescaling
  c_float c;     ///< objective scaling
  c_float cinv;  ///< objective rescaling
} QPALMScaling;


/**
 * Solver return information
 */
typedef struct {
  c_int iter;          ///< number of iterations taken
  c_int iter_out;      ///< number of outer iterations
  char  status[32];    ///< status string, e.g. 'solved'
  c_int status_val;    ///< status as c_int, defined in constants.h

  c_float obj_val;          ///< primal objective
  c_float pri_res_norm;     ///< norm of primal residual
  c_float dua_res_norm;     ///< norm of dual residual
  c_float dua2_res_norm;    ///< norm of intermediate dual residual (minus proximal term)

  #ifdef PROFILING
  c_float setup_time;  ///< time taken for setup phase (seconds)
  c_float solve_time;  ///< time taken for solve phase (seconds)
  c_float run_time;    ///< total time  (seconds)
  #endif

} QPALMInfo;

/**********************************
* Main structures and Data Types *
**********************************/

/**
 * Data structure
 */
typedef struct {
  size_t    n; ///< number of variables n
  size_t   m; ///< number of constraints m
  cholmod_sparse *Q; ///< quadratic part of the cost Q in csc format (size n x n). It
              ///  can be either the full Q or only the upper triangular part. 
  cholmod_sparse *A; ///< linear constraints matrix A in csc format (size m x n)
  c_float *q; ///< dense array for linear part of cost function (size n)
  c_float *bmin; ///< dense array for lower bound (size m)
  c_float *bmax; ///< dense array for upper bound (size m)
} QPALMData;


/**
 * Settings struct
 */
typedef struct {
  c_int   max_iter;      ///< maximum iterations
  c_float eps_abs;       ///< absolute convergence tolerance
  c_float eps_rel;       ///< relative convergence tolerance
  c_float eps_abs_in;    ///< intermediate absolute convergence tolerance
  c_float eps_rel_in;    ///< intermediate relative convergence tolerance
  c_float rho;           ///< tolerance scaling factor 
  c_float eps_prim_inf;  ///< primal infeasibility tolerance
  c_float eps_dual_inf;  ///< dual infeasibility tolerance
  c_float theta;         ///< penalty update criterion parameter
  c_float delta;         ///< penalty update factor
  c_float tau_init;      ///< initial stepsize in backtracking
  c_int   memory;        ///< LBFGS memory
  c_int   proximal;      ///< boolean, proximal method of multipliers
  c_float gamma;         ///< proximal penalty parameter
  c_float gamma_upd;     ///< proximal penalty update factor
  c_float gamma_max;     ///< proximal penalty parameter cap
  c_int   scaling;       ///< scaling iterations, if 0 scaling disabled
  c_int   verbose;       ///< boolean, write out progress
  c_int   warm_start;    ///< boolean, warm start
} QPALMSettings;

/**
 * LBFGS struct
 */
typedef struct {
  c_int    curridx;     ///< Starting position in the buffers
  c_int    currmem;     ///< Number of (valid) elements in buffers
  c_int    reset_lbfgs; ///< boolean, reset lbfgs (after sigma update)
  c_float *s;           ///< s = x-x_prev
  c_float *y;           ///< y = dphi-dphi_prev
  c_float  ys;          ///< y'*s
  c_float *Sbuffer;     ///< Buffer for s vectors
  c_float *Ybuffer;     ///< Buffer for y vectors
  c_float *YSbuffer;    ///< Buffer for ys numbers
  c_float  H;           ///< preconditioning constant in lbfgs
  c_float *alpha;       ///< alpha vector
  c_float *q;           ///< running vector
} QPALMLbfgs;

/**
 * Variables for linear system solving (cholmod)
 */
typedef struct {
  cholmod_common c;
  cholmod_factor *LD;
  cholmod_dense *E_temp;
  cholmod_dense *D_temp;
  cholmod_dense *neg_dphi;
  cholmod_dense *d;
  cholmod_dense *Ad;
  cholmod_dense *Qd;
  cholmod_dense *yh;
  cholmod_dense *Atyh;
  c_int reset_newton;
  c_int *active_constraints;
  c_int *active_constraints_old;
  c_int nb_active_constraints;
  c_int *enter;
  c_int nb_enter;
  c_int *leave;
  c_int nb_leave;
  cholmod_dense *At_scale;
  cholmod_sparse *At_sqrt_sigma;
  cholmod_sparse *Q_AtA;
} QPALMCholmod;

/**
 * QPALM Workspace
 */
typedef struct {
  /// Problem data to work on (possibly scaled)
  QPALMData *data;

  /**
   * @name Iterates
   * @{
   */
  c_float *x;        ///< Iterate x
  c_float *y;        ///< Iterate y
  c_float *Ax;       ///< Scaled A * x
  c_float *Qx;       ///< Scaled Q * x
  c_float *Aty;      ///< Aty (useful for sparing one mat_tpose_vec)
  c_float *x_prev;   ///< Previous x

  /** @} */

  /**
   * @name Workspace variables
   * @{ 
   */
  c_float *temp_m;     ///< Placeholder for vector of size m
  c_float *temp_n;     ///< Placeholder for vector of size n
  c_float *sigma;      ///< Penalty vector
  c_float *Axys;       ///< Ax + y./sigma
  c_float *z;          ///< z vector
  c_float *pri_res;    ///< Primal residual
  c_float *pri_res_in; ///< Intermediate Primal residual
  c_float *yh;         ///< Candidate dual update
  c_float *Atyh;       ///< A'*yh
  c_float *df;         ///< Gradient primal objective (+proximal term)
  c_float *x0;         ///< x0, used in proximal method of multipliers
  c_float *xx0;        ///< x-x0
  c_float *dphi;       ///< Gradient Lagrangian
  c_float *neg_dphi;   ///< -dphi, needed as rhs in Cholmod 
  c_float *dphi_prev;  ///< Previous gradient Lagrangian
  c_float *d;          ///< Step

  /** @} */

  /**
   * @name Linesearch variables
   * @{
   */
  c_float *Qd;        ///< Q*d
  c_float *Ad;        ///< A*d
  c_float *sqrt_sigma; ///< Elementwise sqrt(sigma)
  c_float sqrt_delta; ///< sqrt(penalty update factor)
  c_float eta;        ///< Linesearch parameter
  c_float beta;       ///< Linesearch parameter
  c_float *delta;     ///< Linesearch parameter
  c_float *alpha;     ///< Linesearch parameter
  c_float *temp_2m;   ///< Placeholder for vector of size 2m
  c_float *delta2;    ///< delta.*delta
  c_float *delta_alpha; ///< delta.*alpha
  array_element *s;   ///< alpha./delta
  c_int   *index_L;   ///< Index set L
  c_int   *index_P;   ///< Index set P
  c_int   *index_J;   ///< Index set J

  /** @} */

  /**
   * @name Termination criteria variables
   * @{ 
   */
  c_float eps_pri;    ///< Primal tolerance
  c_float eps_dua;    ///< Dual tolerance
  c_float eps_dua_in; ///< Intermediate dual tolerance
  /** @} */

  /**
   * @name Primal infeasibility variables
   * @{
   */
  c_float *delta_y;   ///< Difference of consecutive dual iterates
  c_float *Atdelta_y; ///< A' * delta_y

  /** @} */

  /**
   * @name Dual infeasibility variables
   * @{
   */
  c_float *delta_x;  ///< Difference of consecutive primal iterates
  c_float *Qdelta_x; ///< Q * delta_x
  c_float *Adelta_x; ///< A * delta_x

  /** @} */

  /**
   * @name Temporary vectors used in scaling
   * @{
   */

  c_float *D_temp;   ///< temporary primal variable scaling vectors
  c_float *E_temp;   ///< temporary constraints scaling vectors


  /** @} */

  QPALMCholmod  *chol;     ///< Cholmod variables
  QPALMLbfgs    *lbfgs;    ///< LBFGS variables
  QPALMSettings *settings; ///< Problem settings
  QPALMScaling  *scaling;  ///< Scaling vectors
  QPALMSolution *solution; ///< Problem solution
  QPALMInfo     *info;     ///< Solver information

  # ifdef PROFILING
  QPALMTimer *timer;       ///< Timer object
  # endif // ifdef PROFILING

} QPALMWorkspace;


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef QPALM_TYPES_H
