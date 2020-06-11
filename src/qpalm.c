/**
 * @file qpalm.c
 * @author Ben Hermans
 * @brief QPALM main solver API.
 * @details This file contains the main functions that can be called by the user.
 * The user can load the default settings, setup the workspace with data and settings,
 * run the solver, and cleanup the workspace afterwards.
 */
# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "qpalm.h"
#include "global_opts.h"
#include "constants.h"
#include "validate.h"
#include "lin_alg.h"
#include "util.h"
#include "scaling.h"
#include "linesearch.h"
#include "termination.h"
#include "solver_interface.h"
#include "newton.h"
#include "nonconvex.h"
#include "iteration.h"

#ifdef USE_LADEL
#include "ladel.h"
#elif defined USE_CHOLMOD
#include "cholmod.h"
#include "cholmod_function.h"
#endif

/**********************
* Main API Functions *
**********************/

void qpalm_set_default_settings(QPALMSettings *settings) {

  settings->max_iter                = MAX_ITER;                /* maximum iterations */
  settings->inner_max_iter          = INNER_MAX_ITER;          /* maximum iterations per subproblem */
  settings->eps_abs                 = (c_float)EPS_ABS;        /* absolute convergence tolerance */
  settings->eps_rel                 = (c_float)EPS_REL;        /* relative convergence tolerance */
  settings->eps_abs_in              = (c_float)EPS_ABS_IN;     /* intermediate absolute convergence tolerance */
  settings->eps_rel_in              = (c_float)EPS_REL_IN;     /* intermediate relative convergence tolerance */
  settings->rho                     = (c_float)RHO;            /* tolerance scaling factor */
  settings->eps_prim_inf            = (c_float)EPS_PRIM_INF;   /* primal infeasibility tolerance */
  settings->eps_dual_inf            = (c_float)EPS_DUAL_INF;   /* dual infeasibility tolerance */
  settings->theta                   = (c_float)THETA;          /* penalty update criterion parameter */
  settings->delta                   = (c_float)DELTA;          /* penalty update factor */
  settings->sigma_max               = (c_float)SIGMA_MAX;      /* penalty parameter cap */
  settings->proximal                = PROXIMAL;                /* boolean, proximal method of multipliers*/
  settings->gamma_init              = (c_float)GAMMA_INIT;     /* proximal penalty parameter */
  settings->gamma_upd               = (c_float)GAMMA_UPD;      /* proximal penalty update factor*/
  settings->gamma_max               = (c_float)GAMMA_MAX;      /* proximal penalty parameter cap*/
  settings->scaling                 = SCALING;                 /* boolean, scaling */
  settings->nonconvex               = NONCONVEX;               /* boolean, nonconvex */
  settings->warm_start              = WARM_START;              /* boolean, warm start solver */
  settings->verbose                 = VERBOSE;                 /* boolean, write out progress */
  settings->print_iter              = PRINT_ITER;              /* frequency of printing */
  settings->reset_newton_iter       = RESET_NEWTON_ITER;       /* frequency of performing a full Cholesky factorization */
  settings->enable_dual_termination = ENABLE_DUAL_TERMINATION; /* allow for dual termination (useful in branch and bound) */
  settings->dual_objective_limit    = DUAL_OBJECTIVE_LIMIT;    /* termination value for the dual objective (useful in branch and bound) */
  settings->time_limit              = TIME_LIMIT;              /* time limit */
  settings->ordering                = ORDERING;                /* ordering */
  settings->factorization_method    = FACTORIZATION_METHOD;    /* factorization method (kkt or schur) */
}


QPALMWorkspace* qpalm_setup(const QPALMData *data, const QPALMSettings *settings) {

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
  work->gamma = work->settings->gamma_init;

  size_t n = data->n;
  size_t m = data->m;

  //Initialize the solver for the linear system
  work->solver = c_calloc(1, sizeof(QPALMSolver));
  solver_common common, *c;
  c = &common;
  #ifdef USE_LADEL
  #elif defined USE_CHOLMOD
  CHOLMOD(start)(c);
  cholmod_set_settings(c);
  #endif

  // Copy problem data into workspace
  work->data       = c_calloc(1, sizeof(QPALMData));
  work->data->n    = data->n;           
  work->data->m    = data->m;                   
  work->data->bmin = vec_copy(data->bmin, m);      
  work->data->bmax = vec_copy(data->bmax, m);       
  work->data->q    = vec_copy(data->q, n);          
  work->data->c    = data->c;

  #ifdef USE_LADEL
  work->data->A = ladel_sparse_allocate_and_copy(data->A); 
  work->data->Q = ladel_sparse_allocate_and_copy(data->Q);
  ladel_to_upper_diag(work->data->Q);
  #elif defined USE_CHOLMOD
  work->data->A    = CHOLMOD(copy_sparse)(data->A, c); 
  work->data->A->stype = 0;   
  work->data->Q    = CHOLMOD(copy_sparse)(data->Q, c);     
  #endif

  qpalm_set_factorization_method(work);

  // Allocate internal solver variables 
  work->x        = c_calloc(n, sizeof(c_float));
  work->y        = c_calloc(m, sizeof(c_float));
  work->Ax       = c_calloc(m, sizeof(c_float));
  work->Qx       = c_calloc(n, sizeof(c_float));
  work->x_prev   = c_calloc(n, sizeof(c_float));
  work->Aty      = c_calloc(n, sizeof(c_float));

  work->x0 = c_calloc(n, sizeof(c_float));

  work->initialized = FALSE;

  // Workspace variables
  work->temp_m   = c_calloc(m, sizeof(c_float));
  work->temp_n   = c_calloc(n, sizeof(c_float));
  work->sigma = c_calloc(m, sizeof(c_float));
  work->sigma_inv = c_calloc(m, sizeof(c_float));
  work->nb_sigma_changed = 0;

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

    #ifdef USE_LADEL
    work->solver->E_temp = c_calloc(m, sizeof(c_float));
    work->E_temp = work->solver->E_temp;
    work->solver->D_temp = c_calloc(n, sizeof(c_float));
    work->D_temp = work->solver->D_temp;
    #elif defined USE_CHOLMOD
   
    // Allocate cholmod_dense pointers to E_temp and D_temp
    work->solver->E_temp = CHOLMOD(allocate_dense)(m, 1, m, CHOLMOD_REAL, c);
    work->E_temp = work->solver->E_temp->x;
    work->solver->D_temp = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, c);
    work->D_temp = work->solver->D_temp->x;
    #endif

    // Scale data
    scale_data(work);
    
  }
  else {
    work->scaling = QPALM_NULL;
  }

  // Solver variables
  work->solver->active_constraints = c_calloc(m, sizeof(c_int));
  work->solver->active_constraints_old = c_calloc(m, sizeof(c_int));
  vec_set_scalar_int(work->solver->active_constraints_old, FALSE, m);
  work->solver->reset_newton = TRUE;
  work->solver->enter = c_calloc(m, sizeof(c_int));
  work->solver->leave = c_calloc(m, sizeof(c_int));

  #ifdef USE_LADEL
  if (work->solver->factorization_method == FACTORIZE_KKT)
  {
    work->solver->rhs_kkt = c_malloc((n+m)*sizeof(c_float));
    work->solver->sol_kkt = c_malloc((n+m)*sizeof(c_float));
    c_int kkt_nzmax = work->data->Q->nzmax + work->data->A->nzmax + m;
    work->solver->kkt_full = ladel_sparse_alloc(n+m, n+m, kkt_nzmax, UPPER, TRUE, FALSE);
    work->solver->kkt = ladel_sparse_alloc(n+m, n+m, kkt_nzmax, UPPER, TRUE, TRUE);
    work->solver->first_row_A = c_malloc(m*sizeof(c_int));
    work->solver->first_elem_A = c_malloc(m*sizeof(c_float));
    work->solver->sym = ladel_symbolics_alloc(m+n);

    c->array_int_ncol1 = work->index_L; /* Avoid allocating full workspace */
    work->solver->At = ladel_transpose(work->data->A, TRUE, c);
    c->array_int_ncol1 = NULL;
  } else if (work->solver->factorization_method == FACTORIZE_SCHUR)
  {

  }
  
  work->solver->neg_dphi = c_calloc(n, sizeof(c_float));
  work->neg_dphi = work->solver->neg_dphi; 
  work->solver->d = c_calloc(n, sizeof(c_float));
  work->d = work->solver->d;
  work->solver->Qd = c_calloc(n, sizeof(c_float));
  work->Qd = work->solver->Qd;
  work->solver->Ad = c_calloc(m, sizeof(c_float));
  work->Ad = work->solver->Ad;
  work->solver->yh = c_calloc(m, sizeof(c_float));
  work->yh = work->solver->yh;
  work->solver->Atyh = c_calloc(n, sizeof(c_float));
  work->Atyh = work->solver->Atyh;
  work->solver->At_scale = c_calloc(m, sizeof(c_float));

  work->solver->first_factorization = TRUE;
  
  if (work->settings->enable_dual_termination)
    work->solver->sym_Q = ladel_symbolics_alloc(n);

  #elif defined USE_CHOLMOD
  work->solver->neg_dphi = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, c);
  work->neg_dphi = work->solver->neg_dphi->x; 
  work->solver->d = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, c);
  work->d = work->solver->d->x;
  work->solver->Qd = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, c);
  work->Qd = work->solver->Qd->x;
  work->solver->Ad = CHOLMOD(allocate_dense)(m, 1, m, CHOLMOD_REAL, c);
  work->Ad = work->solver->Ad->x;
  work->solver->yh = CHOLMOD(allocate_dense)(m, 1, m, CHOLMOD_REAL, c);
  work->yh = work->solver->yh->x;
  work->solver->Atyh = CHOLMOD(allocate_dense)(n, 1, n, CHOLMOD_REAL, c);
  work->Atyh = work->solver->Atyh->x;
  work->solver->At_scale = CHOLMOD(allocate_dense)(m, 1, m, CHOLMOD_REAL, c);

  #endif
  // Set parameters in case the QP is nonconvex
  if (work->settings->nonconvex) {
    set_settings_nonconvex(work, c);
  }

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

  #ifdef USE_LADEL
  #elif defined USE_CHOLMOD
  CHOLMOD(finish)(c);
  #endif

  // Return workspace structure
  return work;
}


void qpalm_warm_start(QPALMWorkspace *work, c_float *x_warm_start, c_float *y_warm_start) {
    
    work->gamma = work->settings->gamma_init;

    // If we have previously solved the problem, then just count the warm start as the setup time
    #ifdef PROFILING
    if (work->info->status_val != QPALM_UNSOLVED) work->info->setup_time = 0; 
    qpalm_tic(work->timer); // Start timer
    #endif /* ifdef PROFILING */

    size_t n = work->data->n;
    size_t m = work->data->m;
    solver_common common, *c;
    c = &common;
    
    #ifdef USE_LADEL
    #elif defined USE_CHOLMOD
    CHOLMOD(start)(c);
    cholmod_set_settings(c);
    #endif

    if (x_warm_start != NULL) {
      prea_vec_copy(x_warm_start, work->x, n);
      // Scale initial vectors x, xprev, x0, if they are warm-started
      if (work->settings->scaling) {
        vec_ew_prod(work->x, work->scaling->Dinv, work->x, n);
      }
      prea_vec_copy(work->x, work->x0, n);
      prea_vec_copy(work->x, work->x_prev, n);
      
      //NB to link to cholmod here neg_dphi and chol->neg_dphi are used as x
      // work->d = work->solver->neg_dphi->x;
      prea_vec_copy(work->x, work->neg_dphi, n);

      mat_vec(work->data->Q, work->solver->neg_dphi, work->solver->Qd, c);
      if (work->settings->proximal) {
        vec_add_scaled(work->Qd, work->x, work->Qx, 1/work->settings->gamma_init, n);
      } else {
        prea_vec_copy(work->Qd, work->Qx, n);
      }
      mat_vec(work->data->A, work->solver->neg_dphi, work->solver->Ad, c);
      prea_vec_copy(work->Ad, work->Ax, m);

      work->info->objective = compute_objective(work);

    } else {
      vec_set_scalar(work->x, 0., n);
      vec_set_scalar(work->x_prev, 0., n);
      vec_set_scalar(work->x0, 0., n);
      vec_set_scalar(work->Qx, 0., n);
      vec_set_scalar(work->Ax, 0., m);
      work->info->objective = 0.0;
    }

    if (y_warm_start != NULL) {
      prea_vec_copy(y_warm_start, work->y, m);
      if (work->settings->scaling) {
        vec_ew_prod(work->y, work->scaling->Einv, work->y, m);
        vec_self_mult_scalar(work->y, work->scaling->c, m);
      }
    } else {
      vec_set_scalar(work->y, 0., m);
    }
    
    initialize_sigma(work, c);

    work->initialized = TRUE;
    
    #ifdef USE_LADEL
    #elif defined USE_CHOLMOD
    CHOLMOD(finish)(c);
    #endif

    #ifdef PROFILING
    work->info->setup_time += qpalm_toc(work->timer); // Start timer
    #endif /* ifdef PROFILING */

}

void qpalm_solve(QPALMWorkspace *work) {
  
  // Print header
  #ifdef PRINTING
  if (work->settings->verbose) {
    print_header();
  }
  #endif

  //Set initial workspace variables
  work->eps_abs_in = work->settings->eps_abs_in;
  work->eps_rel_in = work->settings->eps_rel_in;
  work->solver->reset_newton = TRUE;
  work->gamma = work->settings->gamma_init;
  work->gamma_maxed = (FALSE || work->settings->nonconvex);
  vec_set_scalar_int(work->solver->active_constraints_old, FALSE, work->data->m);

  //Check if the internal variables were correctly initialized. A path that leads to
  //incorrect initialization (and subsequent program crash) is that the user forgets
  //to call qpalm_warm_start while having set settings->warm_start=TRUE, which causes 
  //qpalm_setup to not initialize some variables.
  if (!work->initialized) {
    qpalm_warm_start(work, NULL, NULL);
  }  

  // Start the timer (after warm_start because this is already added to the setup time)
  #ifdef PROFILING
  qpalm_tic(work->timer); // Start timer
  c_float current_time;
  #endif /* ifdef PROFILING */
  
  size_t n = work->data->n;
  size_t m = work->data->m;
  //Initialize CHOLMOD and its settings
  solver_common common, *c;
  c = &common;
  #ifdef USE_LADEL
  if (work->solver->factorization_method == FACTORIZE_KKT)
  {
    c = ladel_workspace_allocate(n+m);
  } else if (work->solver->factorization_method == FACTORIZE_SCHUR)
  {
    c = ladel_workspace_allocate(n);
  }
  ladel_work *c2;
  if (work->settings->enable_dual_termination && work->solver->factorization_method == FACTORIZE_KKT) 
    c2 = ladel_workspace_allocate(n);
  else
    c2 = c;
  #elif defined USE_CHOLMOD
  CHOLMOD(start)(c);
  cholmod_set_settings(c);
  #endif


  //Provide LD factor of Q in case dual_termination is enabled
  //NB use neg_dphi = Aty+q and D_temp = Q^-1(Aty+q) to link with cholmod
  //NB assume Q is positive definite
  if (work->settings->enable_dual_termination) {
    #ifdef USE_LADEL
    if (work->solver->LD_Q) ladel_factor_free(work->solver->LD_Q);
    ladel_factorize(work->data->Q, work->solver->sym_Q, work->settings->ordering, &work->solver->LD_Q, c2);
    work->info->dual_objective = compute_dual_objective(work, c2);    
    #elif defined USE_CHOLMOD
    if (work->solver->LD_Q) CHOLMOD(free_factor)(&work->solver->LD_Q, c);
    work->solver->LD_Q = CHOLMOD(analyze) (work->data->Q, c);
    CHOLMOD(factorize) (work->data->Q, work->solver->LD_Q, c);
    work->info->dual_objective = compute_dual_objective(work, c);    
    #endif
  } else {
    work->info->dual_objective = QPALM_NULL;
  }

  c_int iter;
  c_int iter_out = 0;
  c_int prev_iter = 0; /* iteration number at which the previous subproblem finished*/

  for (iter = 0; iter < work->settings->max_iter; iter++) {

    compute_residuals(work, c);
  
    if (check_termination(work)) {
      work->info->iter = iter;
      work->info->iter_out = iter_out;
      /* Update solve time and run time */
      #ifdef PROFILING
        work->info->solve_time = qpalm_toc(work->timer);
        work->info->run_time = work->info->setup_time +
                           work->info->solve_time;
      #endif /* ifdef PROFILING */
      work->initialized = FALSE;

      #ifdef USE_LADEL
      c = ladel_workspace_free(c);
      if (work->settings->enable_dual_termination) 
        c2 = ladel_workspace_free(c2);
      #elif defined USE_CHOLMOD
      CHOLMOD(finish)(c);
      #endif

      #ifdef PRINTING
      if (work->settings->verbose) {
        work->info->objective = compute_objective(work);
        print_iteration(iter, work); 
        print_final_message(work);
      }
      #endif

      return; 
    } else if (check_subproblem_termination(work)) {
      prea_vec_copy(work->yh, work->y, m);
      prea_vec_copy(work->Atyh, work->Aty, n);

      if(work->settings->enable_dual_termination) {

        work->info->dual_objective = compute_dual_objective(work, c);
        if (work->info->dual_objective > work->settings->dual_objective_limit) {
          
          update_status(work->info, QPALM_DUAL_TERMINATED);
          store_solution(work);

          work->info->iter = iter;
          work->info->iter_out = iter_out;
          /* Update solve time and run time */
          #ifdef PROFILING
            work->info->solve_time = qpalm_toc(work->timer);
            work->info->run_time = work->info->setup_time +
                              work->info->solve_time;
          #endif /* ifdef PROFILING */
          work->initialized = FALSE;

          #ifdef USE_LADEL
          ladel_workspace_free(c);
          if (work->settings->enable_dual_termination) 
            c2 = ladel_workspace_free(c2);
          #elif defined USE_CHOLMOD
          CHOLMOD(finish)(c);
          #endif

          #ifdef PRINTING
          if (work->settings->verbose) {
            work->info->objective = compute_objective(work);
            print_iteration(iter, work); 
            print_final_message(work);
          }
          #endif

          return; 
        }
      }

      work->eps_abs_in = c_max(work->settings->eps_abs, work->settings->rho*work->eps_abs_in);
      work->eps_rel_in = c_max(work->settings->eps_rel, work->settings->rho*work->eps_rel_in);
       
      if (iter_out > 0 && work->info->pri_res_norm > work->eps_pri) {
        update_sigma(work, c);
      } 

      if(work->settings->proximal) {
        if (!work->gamma_maxed && iter_out > 0 && work->solver->nb_enter == 0 && work->solver->nb_leave == 0 && work->info->pri_res_norm < work->eps_pri) {
          //Axys = Ax + y./sigma
            vec_ew_div(work->y, work->sigma, work->temp_m, work->data->m);
            vec_add_scaled(work->Ax, work->temp_m, work->Axys, 1, work->data->m);
            set_active_constraints(work);
            set_entering_leaving_constraints(work);
            if (work->solver->nb_enter == 0 && work->solver->nb_leave == 0) {
              // c_print("Boosting gamma on iter: %d\n", iter);
              boost_gamma(work, c);
            } else {
              update_gamma(work);
            }
          } else {
            update_gamma(work);
          }
        
        prea_vec_copy(work->x, work->x0, work->data->n);
      }

      prea_vec_copy(work->pri_res, work->pri_res_in, m);
      iter_out++;
      prev_iter = iter;

      #ifdef PRINTING
      if (work->settings->verbose && mod(iter, work->settings->print_iter) == 0) {
        c_print("%4ld | ---------------------------------------------------\n", iter);
      }
      #endif
      
    
    } else if (iter == prev_iter + work->settings->inner_max_iter){ 
      
      if (iter_out > 0 && work->info->pri_res_norm > work->eps_pri) {
        update_sigma(work, c);
      } 

      if(work->settings->proximal) {
        update_gamma(work);
        prea_vec_copy(work->x, work->x0, work->data->n);
      }

      prea_vec_copy(work->pri_res, work->pri_res_in, m);
      iter_out++;
      prev_iter = iter;

    } else {

      if (mod(iter, work->settings->reset_newton_iter) == 0) work->solver->reset_newton = TRUE; 
      update_primal_iterate(work, c);

      #ifdef PRINTING
      if (work->settings->verbose && mod(iter, work->settings->print_iter) == 0) {
        work->info->objective = compute_objective(work);
        print_iteration(iter, work);
      }
      #endif

    
    }

    #ifdef PROFILING
    current_time = work->info->setup_time + qpalm_toc(work->timer); // Start timer
    if (current_time > work->settings->time_limit) {
      update_status(work->info, QPALM_TIME_LIMIT_REACHED);
      work->info->iter = iter;
      work->info->iter_out = iter_out;
      store_solution(work);
      #ifdef PROFILING
        work->info->solve_time = qpalm_toc(work->timer);
        work->info->run_time = work->info->setup_time +
                              work->info->solve_time;
      #endif /* ifdef PROFILING */
      work->initialized = FALSE;

      #ifdef USE_LADEL
      ladel_workspace_free(c);
      if (work->settings->enable_dual_termination) 
        c2 = ladel_workspace_free(c2);
      #elif defined USE_CHOLMOD
      CHOLMOD(finish)(c);
      #endif

      #ifdef PRINTING
        if (work->settings->verbose)
          print_final_message(work);
      #endif /* PRINTING */

      return;
    }

    #endif /* ifdef PROFILING */
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
  work->initialized = FALSE;

  #ifdef USE_LADEL
  ladel_workspace_free(c);
  if (work->settings->enable_dual_termination) 
        c2 = ladel_workspace_free(c2);
  #elif defined USE_CHOLMOD
  CHOLMOD(finish)(c);
  #endif

  #ifdef PRINTING
    if (work->settings->verbose)
      print_final_message(work);
  #endif /* PRINTING */
}


void qpalm_update_settings(QPALMWorkspace* work, const QPALMSettings *settings) {
  // Validate settings
  if (!validate_settings(settings)) {
    # ifdef PRINTING
      c_eprint("Settings validation returned failure");
    # endif /* ifdef PRINTING */
    update_status(work->info, QPALM_ERROR);
    return;
  }
  if (work->settings->scaling > settings->scaling) {
    # ifdef PRINTING
      c_eprint("Decreasing the number of scaling iterations is not allowed");
    # endif /* ifdef PRINTING */
    update_status(work->info, QPALM_ERROR);
    return;
  } else if (work->settings->scaling < settings->scaling) {
    // Save current scaling vectors
    prea_vec_copy(work->scaling->D, work->temp_n, work->data->n);
    prea_vec_copy(work->scaling->E, work->temp_m, work->data->m);
    c_float c_temp = work->scaling->c;

    // Perform the remaining scaling iterations
    work->settings->scaling = settings->scaling - work->settings->scaling;
    scale_data(work);
    
    // Compute the total scaling vectors
    vec_ew_prod(work->scaling->D, work->temp_n, work->scaling->D, work->data->n);
    vec_ew_prod(work->scaling->E, work->temp_m, work->scaling->E, work->data->m);
    work->scaling->c *= c_temp;
    // Save the inverses
    vec_ew_recipr(work->scaling->D, work->scaling->Dinv, work->data->n);
    vec_ew_recipr(work->scaling->E, work->scaling->Einv, work->data->m);
    work->scaling->cinv = 1/work->scaling->c;

    #ifdef USE_LADEL
    work->solver->first_factorization = TRUE;
    solver_common common, *c;
    c = &common;
    if (work->solver->factorization_method == FACTORIZE_KKT)
    {
      c->array_int_ncol1 = work->index_L; 
      work->solver->At = ladel_sparse_free(work->solver->At);
      work->solver->At = ladel_transpose(work->data->A, TRUE, c);
      c->array_int_ncol1 = NULL;
    }
    #endif
  }
 
  // Copy settings
  c_free(work->settings);
  work->settings = copy_settings(settings);
  work->sqrt_delta = c_sqrt(work->settings->delta);
}

void qpalm_update_bounds(QPALMWorkspace *work, const c_float *bmin, const c_float *bmax) {
  // Validate bounds
  size_t j;
  size_t m = work->data->m;
  if (bmin != NULL && bmax != NULL) {
    for (j = 0; j < m; j++) {
      if (bmin[j] > bmax[j]) {
        # ifdef PRINTING
          c_eprint("Lower bound at index %d is greater than upper bound: %.4e > %.4e",
                  (int)j, work->data->bmin[j], work->data->bmax[j]);
        # endif /* ifdef PRINTING */
        update_status(work->info, QPALM_ERROR);
        return;
      }
    }
  }
  
  if (bmin != NULL) {
      prea_vec_copy(bmin, work->data->bmin, m); 
  }
  if (bmax != NULL) {
      prea_vec_copy(bmax, work->data->bmax, m);
  }     
  
  if (work->settings->scaling) {
    if (bmin != NULL) {
      vec_ew_prod(work->scaling->E, work->data->bmin, work->data->bmin, m);
    }
    if (bmax != NULL) {
      vec_ew_prod(work->scaling->E, work->data->bmax, work->data->bmax, m);
    } 
  }
}

void qpalm_update_q(QPALMWorkspace *work, const c_float *q) {
  size_t n = work->data->n;
  prea_vec_copy(q, work->data->q, n);    
  if (work->settings->scaling) {
    vec_ew_prod(work->scaling->D, work->data->q, work->data->q, n);    
    // Update cost scaling scalar
    c_float c_old = work->scaling->c, c_ratio;
    if (work->settings->proximal) {
      vec_add_scaled(work->Qx, work->x, work->Qx, -1/work->gamma, work->data->n);
    }
    vec_add_scaled(work->data->q, work->Qx, work->temp_n, work->scaling->cinv, n);
    work->scaling->c = 1/c_max(1.0, vec_norm_inf(work->temp_n, n));
    work->scaling->cinv = 1/work->scaling->c;
    vec_self_mult_scalar(work->data->q, work->scaling->c, n);

    c_ratio = work->scaling->c/c_old;

    solver_common common, *c;
    c = &common;
    #ifdef USE_LADEL
    ladel_scale_scalar(work->data->Q, c_ratio);
    if (work->solver->factorization_method == FACTORIZE_KKT)
    {
      ladel_int index;
      for (index = 0; index < work->data->Q->nzmax; index++) 
        work->solver->kkt->x[index] *= c_ratio;
    }
    work->solver->reset_newton = TRUE;
    #elif defined USE_CHOLMOD
    CHOLMOD(start)(c);
    cholmod_dense *scalar = CHOLMOD(ones)(1,1,CHOLMOD_REAL, c);
    c_float *scalarx = scalar->x;
    scalarx[0] = work->scaling->c/c_old;
    CHOLMOD(scale)(scalar, CHOLMOD_SCALAR, work->data->Q, c);
    CHOLMOD(free_dense)(&scalar, c);
    CHOLMOD(finish)(c);
    #endif

    vec_self_mult_scalar(work->Qx, c_ratio, n);
    if (work->settings->proximal) {
      work->gamma = work->settings->gamma_init;
      vec_add_scaled(work->Qx, work->x, work->Qx, 1/work->gamma, work->data->n);    
    }
  }
}


void qpalm_cleanup(QPALMWorkspace *work) {
  
  if (work) { // If workspace has been allocated
    // Free Data
    solver_common common, *c;
    c = &common;
    if (work->data) {
      #ifdef USE_LADEL
      work->data->Q = ladel_sparse_free(work->data->Q);

      work->data->A = ladel_sparse_free(work->data->A);
      #elif defined USE_CHOLMOD
      CHOLMOD(start)(c);

      if (work->data->Q) CHOLMOD(free_sparse)(&work->data->Q, c);

      if (work->data->A) CHOLMOD(free_sparse)(&work->data->A, c);

      CHOLMOD(finish)(c);
      #endif

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
    
    if (work->sigma_inv) c_free(work->sigma_inv);

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
    if (work->solver) {
      if (work->solver->active_constraints) c_free(work->solver->active_constraints);

      if (work->solver->active_constraints_old) c_free(work->solver->active_constraints_old);

      if (work->solver->enter) c_free(work->solver->enter);

      if (work->solver->leave) c_free(work->solver->leave);

      #ifdef USE_LADEL

      work->solver->sol_kkt = ladel_free(work->solver->sol_kkt);

      work->solver->rhs_kkt = ladel_free(work->solver->rhs_kkt);

      work->solver->D_temp = ladel_free(work->solver->D_temp);

      work->solver->E_temp = ladel_free(work->solver->E_temp);

      work->solver->neg_dphi = ladel_free(work->solver->neg_dphi);

      work->solver->d = ladel_free(work->solver->d);

      work->solver->Qd = ladel_free(work->solver->Qd);

      work->solver->Ad = ladel_free(work->solver->Ad);

      work->solver->yh = ladel_free(work->solver->yh);

      work->solver->Atyh = ladel_free(work->solver->Atyh);

      work->solver->LD = ladel_factor_free(work->solver->LD);

      work->solver->LD_Q = ladel_factor_free(work->solver->LD_Q);

      work->solver->sym = ladel_symbolics_free(work->solver->sym);

      work->solver->sym_Q = ladel_symbolics_free(work->solver->sym_Q);

      work->solver->At_scale = ladel_free(work->solver->At_scale);

      work->solver->At_sqrt_sigma = ladel_sparse_free(work->solver->At_sqrt_sigma);

      work->solver->At = ladel_sparse_free(work->solver->At);

      work->solver->kkt = ladel_sparse_free(work->solver->kkt);

      work->solver->kkt_full = ladel_sparse_free(work->solver->kkt_full);

      work->solver->first_row_A = ladel_free(work->solver->first_row_A);

      work->solver->first_elem_A = ladel_free(work->solver->first_elem_A);

      #elif defined USE_CHOLMOD

      CHOLMOD(start)(c);

      if (work->solver->D_temp) CHOLMOD(free_dense)(&work->solver->D_temp, c);

      if (work->solver->E_temp) CHOLMOD(free_dense)(&work->solver->E_temp, c);

      if (work->solver->neg_dphi) CHOLMOD(free_dense)(&work->solver->neg_dphi, c);

      if (work->solver->d) CHOLMOD(free_dense)(&work->solver->d, c);

      if (work->solver->Qd) CHOLMOD(free_dense)(&work->solver->Qd, c);

      if (work->solver->Ad) CHOLMOD(free_dense)(&work->solver->Ad, c);

      if (work->solver->yh) CHOLMOD(free_dense)(&work->solver->yh, c);

      if (work->solver->Atyh) CHOLMOD(free_dense)(&work->solver->Atyh, c);

      if (work->solver->LD) CHOLMOD(free_factor)(&work->solver->LD, c);

      if (work->solver->LD_Q) CHOLMOD(free_factor)(&work->solver->LD_Q, c);

      if (work->solver->At_scale) CHOLMOD(free_dense)(&work->solver->At_scale, c);

      if (work->solver->At_sqrt_sigma) CHOLMOD(free_sparse)(&work->solver->At_sqrt_sigma, c);

      CHOLMOD(finish)(c);
      #endif

      c_free(work->solver);      
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


# ifdef __cplusplus
}
# endif // ifdef __cplusplus