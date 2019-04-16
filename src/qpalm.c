/**
 * @file qpalm.c
 * @author Ben Hermans
 * @brief QPALM main solver API.
 * @details This file contains the main functions that can be called by the user.
 * The user can load the default settings, setup the workspace with data and settings,
 * run the solver, and cleanup the workspace afterwards.
 */

#include "qpalm.h"
#include "constants.h"
#include "validate.h"
#include "lin_alg.h"
#include "util.h"
#include "scaling.h"
#include "linesearch.h"
#include "termination.h"
#include "cholmod.h"
#include "cholmod_function.h"
#include "cholmod_interface.h"
#include "newton.h"

/**********************
* Main API Functions *
**********************/
void qpalm_set_default_settings(QPALMSettings *settings) {

  settings->max_iter      = MAX_ITER;                /* maximum iterations to take */
  settings->eps_abs       = (c_float)EPS_ABS;        /* absolute convergence tolerance */
  settings->eps_rel       = (c_float)EPS_REL;        /* relative convergence tolerance */
  settings->eps_abs_in    = (c_float)EPS_ABS_IN;     /* intermediate absolute convergence tolerance */
  settings->eps_rel_in    = (c_float)EPS_REL_IN;     /* intermediate relative convergence tolerance */
  settings->rho           = (c_float)RHO;            /* tolerance scaling factor */
  settings->eps_prim_inf  = (c_float)EPS_PRIM_INF;   /* primal infeasibility tolerance */
  settings->eps_dual_inf  = (c_float)EPS_DUAL_INF;   /* dual infeasibility tolerance */
  settings->theta         = (c_float)THETA;          /* penalty update criterion parameter */
  settings->delta         = (c_float)DELTA;          /* penalty update factor */
  settings->tau_init      = (c_float)TAU_INIT;       /* initial stepsize for backtracking */
  settings->proximal      = PROXIMAL;                /* boolean, proximal method of multipliers*/
  settings->gamma         = (c_float)GAMMA;          /* proximal penalty parameter */
  settings->gamma_upd     = (c_float)GAMMA_UPD;      /* proximal penalty update factor*/
  settings->gamma_max     = (c_float)GAMMA_MAX;      /* proximal penalty parameter cap*/
  settings->scaling       = SCALING;                 /* boolean, scaling */
  settings->warm_start    = WARM_START;              /* boolean, warm start solver */
  settings->verbose       = VERBOSE;                 /* boolean, write out progress */
}


QPALMWorkspace* qpalm_setup(const QPALMData *data, QPALMSettings *settings, cholmod_common *c) {
  QPALMWorkspace *work; // Workspace

  // Validate data
  if (!validate_data(data)) {
# ifdef PRINTING
    c_eprint("Data validation returned failure");
# endif /* ifdef PRINTING */
    return QPALM_NULL;
  }

  // Validate settings
  if (!validate_settings(settings)) {
# ifdef PRINTING
    c_eprint("Settings validation returned failure");
# endif /* ifdef PRINTING */
    return QPALM_NULL;
  }

  // Allocate empty workspace
  work = c_calloc(1, sizeof(QPALMWorkspace));

  if (!work) {
# ifdef PRINTING
    c_eprint("allocating work failure");
# endif /* ifdef PRINTING */
    return QPALM_NULL;
  }

  // Start and allocate directly timer
  # ifdef PROFILING
    work->timer = c_malloc(sizeof(QPALMTimer));
    qpalm_tic(work->timer);
  # endif /* ifdef PROFILING */

  // Copy settings
  work->settings = copy_settings(settings);
  work->sqrt_delta = c_sqrt(work->settings->delta);

  size_t n = data->n;
  size_t m = data->m;

  //Initialize CHOLMOD and its settings
  work->chol = c_calloc(1, sizeof(QPALMCholmod));
  work->chol->c = *c;
  CHOLMOD(start)(&work->chol->c);
  cholmod_set_settings(&work->chol->c);

  // Copy problem data into workspace
  work->data       = c_calloc(1, sizeof(QPALMData));
  work->data->n    = data->n;                             // Number of variables
  work->data->m    = data->m;                             // Number of linear constraints
  work->data->bmin = vec_copy(data->bmin, m);       // Lower bounds on constraints
  work->data->bmax = vec_copy(data->bmax, m);       // Upper bounds on constraints
  work->data->q    = vec_copy(data->q, n);          // Linear part of cost function

  work->data->A    = CHOLMOD(copy_sparse)(data->A, &work->chol->c);               // Linear constraints matrix
  work->data->Q    = CHOLMOD(copy_sparse)(data->Q, &work->chol->c);                // Cost function matrix

  // Allocate internal solver variables 
  work->x        = c_calloc(n, sizeof(c_float));
  work->y        = c_calloc(m, sizeof(c_float));
  work->Ax       = c_calloc(m, sizeof(c_float));
  work->Qx       = c_calloc(n, sizeof(c_float));
  work->x_prev   = c_calloc(n, sizeof(c_float));
  work->Aty      = c_calloc(n, sizeof(c_float));

  work->x0 = c_calloc(n, sizeof(c_float));

  // Workspace variables
  work->temp_m   = c_calloc(m, sizeof(c_float));
  work->temp_n   = c_calloc(n, sizeof(c_float));
  work->sigma = c_calloc(m, sizeof(c_float));
    
  // If warm-start will not be used, cold start variables x, x0, x_prev, y, Qx, Ax and init sigma
  if (!work->settings->warm_start) {
    warm_start(work, NULL, NULL);
  }

  work->z  = c_calloc(m, sizeof(c_float));
  work->Axys = c_calloc(m, sizeof(c_float));
  work->pri_res = c_calloc(m, sizeof(c_float));
  work->pri_res_in = c_calloc(m, sizeof(c_float));
  work->df = c_calloc(n, sizeof(c_float));
  
  work->xx0 = c_calloc(n, sizeof(c_float));
  work->dphi = c_calloc(n, sizeof(c_float));
  work->dphi_prev = c_calloc(n, sizeof(c_float));

  // Linesearch variables
  work->sqrt_sigma  = c_calloc(m, sizeof(c_float));
  work->delta       = c_calloc(m*2, sizeof(c_float));
  work->alpha       = c_calloc(m*2, sizeof(c_float));
  work->delta2      = c_calloc(m*2, sizeof(c_float));
  work->delta_alpha = c_calloc(m*2, sizeof(c_float));
  work->temp_2m     = c_calloc(m*2, sizeof(c_float));
  work->s           = c_calloc(m*2, sizeof(array_element));
  work->index_L     = c_calloc(m*2, sizeof(c_int));
  work->index_P     = c_calloc(m*2, sizeof(c_int));
  work->index_J     = c_calloc(m*2, sizeof(c_int));

  // Primal infeasibility variables
  work->delta_y   = c_calloc(m, sizeof(c_float));
  work->Atdelta_y = c_calloc(n, sizeof(c_float));

  // Dual infeasibility variables
  work->delta_x  = c_calloc(n, sizeof(c_float));
  work->Qdelta_x = c_calloc(n, sizeof(c_float));
  work->Adelta_x = c_calloc(m, sizeof(c_float));

  // Perform scaling
  if (settings->scaling) {
    // Allocate scaling structure
    work->scaling       = c_malloc(sizeof(QPALMScaling));
    work->scaling->D    = c_calloc(n, sizeof(c_float));
    work->scaling->Dinv = c_calloc(n, sizeof(c_float));
    work->scaling->E    = c_calloc(m, sizeof(c_float));
    work->scaling->Einv = c_calloc(m, sizeof(c_float));

    // Allocate cholmod_dense pointers to E_temp and D_temp
    work->chol->E_temp = CHOLMOD(allocate_dense)(m, 1, m, CHOLMOD_REAL, &work->chol->c);
    work->E_temp = work->chol->E_temp->x;
    work->chol->D_temp = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, &work->chol->c);
    work->D_temp = work->chol->D_temp->x;

    // Scale data
    scale_data(work);
    
  }
  else {
    work->scaling = QPALM_NULL;
  }

  // CHOLMOD variables
  work->chol->neg_dphi = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, &work->chol->c);
  work->neg_dphi = work->chol->neg_dphi->x; 
  work->chol->Qd = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, &work->chol->c);
  work->Qd = work->chol->Qd->x;
  work->chol->Ad = CHOLMOD(allocate_dense)(m, 1, m, CHOLMOD_REAL, &work->chol->c);
  work->Ad = work->chol->Ad->x;
  work->chol->yh = CHOLMOD(allocate_dense)(m, 1, m, CHOLMOD_REAL, &work->chol->c);
  work->yh = work->chol->yh->x;
  work->chol->Atyh = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, &work->chol->c);
  work->Atyh = work->chol->Atyh->x;
  work->chol->active_constraints = c_calloc(m, sizeof(c_int));
  work->chol->active_constraints_old = c_calloc(m, sizeof(c_int));
  vec_set_scalar_int(work->chol->active_constraints_old, FALSE, m);
  work->chol->reset_newton = TRUE;
  work->chol->enter = c_calloc(m, sizeof(c_int));
  work->chol->leave = c_calloc(m, sizeof(c_int));
  work->chol->At_scale = CHOLMOD(allocate_dense)(m, 1, m, CHOLMOD_REAL, &work->chol->c);
  vec_ew_sqrt(work->sigma, work->sqrt_sigma, work->data->m);
  prea_vec_copy(work->sqrt_sigma, work->chol->At_scale->x, work->data->m);
  work->chol->At_sqrt_sigma = CHOLMOD(transpose)(work->data->A, 1, &work->chol->c);
  CHOLMOD(scale)(work->chol->At_scale, CHOLMOD_COL, work->chol->At_sqrt_sigma, &work->chol->c);

  // Allocate solution
  work->solution    = c_calloc(1, sizeof(QPALMSolution));
  work->solution->x = c_calloc(1, n * sizeof(c_float));
  work->solution->y = c_calloc(1, m * sizeof(c_float));
  // Allocate and initialize information
  work->info                = c_calloc(1, sizeof(QPALMInfo));
  update_status(work->info, QPALM_UNSOLVED);
  # ifdef PROFILING
  work->info->solve_time  = 0.0;                    // Solve time to zero
  work->info->run_time    = 0.0;                    // Total run time to zero
  work->info->setup_time  = qpalm_toc(work->timer); // Update timer information
  # endif /* ifdef PROFILING */

  //Finish cholmod
  CHOLMOD(finish)(&work->chol->c);

  // Return workspace structure
  return work;
}


void qpalm_solve(QPALMWorkspace *work) {

  #ifdef PROFILING
  qpalm_tic(work->timer); // Start timer
  #endif /* ifdef PROFILING */
  
  //Initialize CHOLMOD and its settings
  CHOLMOD(start)(&work->chol->c);
  cholmod_set_settings(&work->chol->c);
  
  size_t n = work->data->n;
  size_t m = work->data->m;

  c_int iter;
  c_int iter_out = 0;
  c_float tau;

  for (iter = 0; iter < work->settings->max_iter; iter++) {

    //Axys = Ax + y./sigma
    vec_ew_div(work->y, work->sigma, work->temp_m, m);
    vec_add_scaled(work->Ax, work->temp_m, work->Axys, 1, m);
    //z = min(max(Axys,bmin),bmax)
    vec_ew_mid_vec(work->Axys, work->data->bmin, work->data->bmax, work->z, m);
    //pri_res = Ax-z
    vec_add_scaled(work->Ax, work->z, work->pri_res, -1, m);
    //yh = y + pri_res.*sigma
    vec_ew_prod(work->pri_res, work->sigma, work->temp_m, m);
    vec_add_scaled(work->y, work->temp_m, work->yh, 1, m);
    //df = Qx + q
    vec_add_scaled(work->Qx, work->data->q, work->df, 1, n);
    if (work->settings->proximal) {
      //df = Qx + q +1/gamma*(x-x0)
      // NB work->Qx contains Qx+1/gamma*x
      vec_add_scaled(work->df, work->x0, work->df, -1/work->settings->gamma, n);
    }
    // Atyh = A'*yh
    mat_tpose_vec(work->data->A, work->chol->yh, work->chol->Atyh, &work->chol->c);
    //dphi = df+Atyh
    vec_add_scaled(work->df, work->Atyh, work->dphi, 1, n);
  
    if (check_termination(work)) {
      work->info->iter = iter;
      work->info->iter_out = iter_out;
      /* Update solve time and run time */
      #ifdef PROFILING
        work->info->solve_time = qpalm_toc(work->timer);
        work->info->run_time = work->info->setup_time +
                           work->info->solve_time;
      #endif /* ifdef PROFILING */
      CHOLMOD(finish)(&work->chol->c);

      return; 
    } else if (check_subproblem_termination(work)) {
      prea_vec_copy(work->yh, work->y, m);
      prea_vec_copy(work->Atyh, work->Aty, n);
      work->settings->eps_abs_in = c_max(work->settings->eps_abs,
                                        work->settings->rho*work->settings->eps_abs_in);
      work->settings->eps_rel_in = c_max(work->settings->eps_rel,
                                        work->settings->rho*work->settings->eps_rel_in);
      if (iter_out > 0 && work->info->pri_res_norm > work->eps_pri) {
        c_float *At_scalex = work->chol->At_scale->x;
        for (size_t k = 0; k < m; k++) {
          if (c_absval(work->pri_res[k]) > work->settings->theta*c_absval(work->pri_res_in[k])) {
            work->sigma[k] *= work->settings->delta;
            work->sqrt_sigma[k] *= work->sqrt_delta;
            At_scalex[k] = work->sqrt_delta;
          } else {
            At_scalex[k] = 1.0;
          }
        }
        CHOLMOD(scale)(work->chol->At_scale, CHOLMOD_COL, work->chol->At_sqrt_sigma, &work->chol->c);
      }
      prea_vec_copy(work->pri_res, work->pri_res_in, m);
      
      work->chol->reset_newton = TRUE;
      vec_set_scalar_int(work->chol->active_constraints_old, FALSE, work->data->m);

      if(work->settings->proximal) {
        work->settings->gamma = c_min(work->settings->gamma*work->settings->gamma_upd, work->settings->gamma_max);
        prea_vec_copy(work->x, work->x0, n);
      }

      iter_out++;
    } else {
      // d=newton()
      newton_set_direction(work);
      
      tau = exact_linesearch(work);

      //x_prev = x
      prea_vec_copy(work->x, work->x_prev, n);
      //dphi_prev = dphi 
      prea_vec_copy(work->dphi, work->dphi_prev, n);
      //x = x+tau*d
      vec_add_scaled(work->x, work->d, work->x, tau, n);
      vec_mult_scalar(work->Qd, tau, n); //Qdx used in dua_infeas check
      vec_mult_scalar(work->Ad, tau, m); //Adx used in dua_infeas check
      vec_add_scaled(work->Qx, work->Qd, work->Qx, 1, n);
      vec_add_scaled(work->Ax, work->Ad, work->Ax, 1, m);
    }
  }

  // maxiter reached
  update_status(work->info, QPALM_MAX_ITER_REACHED);
  work->info->iter = iter;
  work->info->iter_out = iter_out;
  store_solution(work);
  #ifdef PROFILING
    work->info->solve_time = qpalm_toc(work->timer);
    work->info->run_time = work->info->setup_time +
                           work->info->solve_time;
  #endif /* ifdef PROFILING */

  CHOLMOD(finish)(&work->chol->c);
}


void qpalm_cleanup(QPALMWorkspace *work) {
  
  if (work) { // If workspace has been allocated
    // Free Data
    if (work->data) {
      CHOLMOD(start)(&work->chol->c);

      if (work->data->Q) CHOLMOD(free_sparse)(&work->data->Q, &work->chol->c);

      if (work->data->A) CHOLMOD(free_sparse)(&work->data->A, &work->chol->c);

      CHOLMOD(finish)(&work->chol->c);

      if (work->data->q) c_free(work->data->q);

      if (work->data->bmin) c_free(work->data->bmin);

      if (work->data->bmax) c_free(work->data->bmax);
      c_free(work->data);
    }

    // Free scaling
    if (work->settings->scaling) {
      if (work->scaling->D) c_free(work->scaling->D);

      if (work->scaling->Dinv) c_free(work->scaling->Dinv);

      if (work->scaling->E) c_free(work->scaling->E);

      if (work->scaling->Einv) c_free(work->scaling->Einv);

      c_free(work->scaling);
    }

    // Free other Variables
    if (work->x) c_free(work->x);
    
    if (work->y) c_free(work->y);

    if (work->Ax) c_free(work->Ax);

    if (work->Qx) c_free(work->Qx);
    
    if (work->x_prev) c_free(work->x_prev);

    if (work->Aty) c_free(work->Aty);

    if (work->temp_m) c_free(work->temp_m);

    if (work->temp_n) c_free(work->temp_n);

    if (work->sigma) c_free(work->sigma);

    if (work->z) c_free(work->z);

    if (work->Axys) c_free(work->Axys);

    if (work->pri_res) c_free(work->pri_res);

    if (work->pri_res_in) c_free(work->pri_res_in);

    if (work->df) c_free(work->df);

    if (work->x0) c_free(work->x0);

    if (work->xx0) c_free(work->xx0);

    if (work->dphi) c_free(work->dphi);

    if (work->dphi_prev) c_free(work->dphi_prev);

    if (work->sqrt_sigma) c_free(work->sqrt_sigma);

    if (work->delta) c_free(work->delta);

    if (work->alpha) c_free(work->alpha);

    if (work->delta2) c_free(work->delta2);

    if (work->delta_alpha) c_free(work->delta_alpha);

    if (work->temp_2m) c_free(work->temp_2m);

    if (work->s) c_free(work->s);

    if (work->index_L) c_free(work->index_L);

    if (work->index_P) c_free(work->index_P);

    if (work->index_J) c_free(work->index_J);

    if (work->delta_y) c_free(work->delta_y);

    if (work->Atdelta_y) c_free(work->Atdelta_y);

    if (work->delta_x) c_free(work->delta_x);

    if (work->Qdelta_x) c_free(work->Qdelta_x);

    if (work->Adelta_x) c_free(work->Adelta_x);

    // Free Settings
    if (work->settings) c_free(work->settings);

    //Free chol struct
    if (work->chol) {
      CHOLMOD(start)(&work->chol->c);

      if (work->chol->D_temp) CHOLMOD(free_dense)(&work->chol->D_temp, &work->chol->c);

      if (work->chol->E_temp) CHOLMOD(free_dense)(&work->chol->E_temp, &work->chol->c);

      if (work->chol->neg_dphi) CHOLMOD(free_dense)(&work->chol->neg_dphi, &work->chol->c);

      if (work->chol->d) CHOLMOD(free_dense)(&work->chol->d, &work->chol->c);

      if (work->chol->Qd) CHOLMOD(free_dense)(&work->chol->Qd, &work->chol->c);

      if (work->chol->Ad) CHOLMOD(free_dense)(&work->chol->Ad, &work->chol->c);

      if (work->chol->yh) CHOLMOD(free_dense)(&work->chol->yh, &work->chol->c);

      if (work->chol->Atyh) CHOLMOD(free_dense)(&work->chol->Atyh, &work->chol->c);

      if (work->chol->LD) CHOLMOD(free_factor)(&work->chol->LD, &work->chol->c);

      if (work->chol->active_constraints) c_free(work->chol->active_constraints);

      if (work->chol->active_constraints_old) c_free(work->chol->active_constraints_old);

      if (work->chol->enter) c_free(work->chol->enter);

      if (work->chol->leave) c_free(work->chol->leave);

      if (work->chol->At_scale) CHOLMOD(free_dense)(&work->chol->At_scale, &work->chol->c);

      if (work->chol->At_sqrt_sigma) CHOLMOD(free_sparse)(&work->chol->At_sqrt_sigma, &work->chol->c);

      CHOLMOD(finish)(&work->chol->c);

      c_free(work->chol);      
    }
    
    // Free solution
    if (work->solution) {
      if (work->solution->x) c_free(work->solution->x);

      if (work->solution->y) c_free(work->solution->y);
      c_free(work->solution);
    }

    // Free timer
    # ifdef PROFILING
      if (work->timer) c_free(work->timer);
    # endif /* ifdef PROFILING */

    // Free information
    if (work->info) c_free(work->info);

    // Free work
    c_free(work);
  }

}