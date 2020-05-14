/**
 * @file cholmod_interface.c
 * @author Ben Hermans
 * @brief Interface and wrapper to cholmod functions
 * @details This file includes all calls to cholmod functions apart from scaling in scaling.c and memory
 * allocation/deallocation in the main functions in qpalm.c. It includes all matrix operations, such as
 * matrix vector products, row- and columnwise norms, cholesky factorizations, factorization updates and
 * solving the linear system. Finally, all the settings relevant to cholmod (and suitesparse) are included
 * in this file as well.
 */

#include "lin_alg.h"
#include "cholmod_interface.h"
#include <stdio.h>

#ifdef USE_LADEL
void mat_vec(solver_sparse *A, solver_dense *x, solver_dense *y, solver_common *c) 
{

}

void mat_tpose_vec(solver_sparse *A, solver_dense *x, solver_dense *y, solver_common *c)
{

}

void mat_inf_norm_cols(solver_sparse *M, c_float *E)
{

}

void mat_inf_norm_rows(solver_sparse *M, c_float *E) 
{

}

#elif defined USE_CHOLMOD

void mat_vec(solver_sparse *A, solver_dense *x, solver_dense *y, solver_common *c) {
    double one [2] = {1,0};
    double zero [2] = {0,0};
    if (x!=y) {
      CHOLMOD(sdmult)(A, 0, one, zero, x, y, c);
    } else {
      cholmod_dense* x2 = CHOLMOD(copy_dense)(x, c);
      CHOLMOD(sdmult)(A, 0, one, zero, x2, y, c);
      CHOLMOD(free_dense)(&x2, c);
    }
}

void mat_tpose_vec(solver_sparse *A, solver_dense *x, solver_dense *y, solver_common *c) {
    double one [2] = {1,0};
    double zero [2] = {0,0};
    if (x!=y) {
      CHOLMOD(sdmult)(A, 1, one, zero, x, y, c);
    } else {
      cholmod_dense* x2 = CHOLMOD(copy_dense)(x, c);
      CHOLMOD(sdmult)(A, 1, one, zero, x2, y, c);
      CHOLMOD(free_dense)(&x2, c);
    }
}

void mat_inf_norm_cols(solver_sparse *M, c_float *E) {
  size_t j;
  c_int k;
  c_float *Mx = M->x;
  c_int *Mp = M->p;

  // Initialize zero max elements
  for (j = 0; j < M->ncol; j++) {
    E[j] = 0.;
  }

  // Compute maximum across columns
  for (j = 0; j < M->ncol; j++) {
    for (k = Mp[j]; k < Mp[j + 1]; k++) {
      E[j] = c_max(c_absval(Mx[k]), E[j]);
    }
  }
}

void mat_inf_norm_rows(solver_sparse *M, c_float *E) {
  size_t j;
  c_int i, k;
  c_int *Mp = M->p;
  c_int *Mi = M->i;
  c_float *Mx = M->x;
 
  // Initialize zero max elements
  for (j = 0; j < M->nrow; j++) {
    E[j] = 0.;
  }

  // Compute maximum across rows
  for (j = 0; j < M->ncol; j++) {
    for (k = Mp[j]; k < Mp[j + 1]; k++) {
      i    = Mi[k];
      E[i] = c_max(c_absval(Mx[k]), E[i]);
    }
  }
}

void ldlchol(cholmod_sparse *M, QPALMWorkspace *work, solver_common *c) {
  if (work->solver->LD) {
      CHOLMOD(free_factor)(&work->solver->LD, c);
  }
  work->solver->LD = CHOLMOD(analyze) (M, c) ;
  if (work->settings->proximal) {
    double beta [2] = {1.0/work->gamma,0};
    CHOLMOD(factorize_p) (M, beta, NULL, 0, work->solver->LD, c);
  } else {
    CHOLMOD(factorize) (M, work->solver->LD, c);
  }
  /* If integers are used, supernodal might fail, so check for this and switch to simplicial if necessary. */
  #ifndef DLONG
  if ((c)->status != CHOLMOD_OK) {
      (c)->supernodal = CHOLMOD_SIMPLICIAL;
      if (work->settings->proximal) {
        double beta [2] = {1.0/work->gamma,0};
        CHOLMOD(factorize_p) (M, beta, NULL, 0, work->solver->LD, c);
      } else {
        CHOLMOD(factorize) (M, work->solver->LD, c);
      }
  }
  #endif /* DLONG */

}

void ldlcholQAtsigmaA(QPALMWorkspace *work, solver_common *c) {
  cholmod_sparse *AtsigmaA;
  cholmod_sparse *QAtsigmaA;
  size_t nb_active = 0;
  for (size_t i = 0; i < work->data->m; i++) {
      if (work->solver->active_constraints[i]){
          work->solver->enter[nb_active] = (c_int)i;
          nb_active++;
      }      
  }
  AtsigmaA = CHOLMOD(aat)(work->solver->At_sqrt_sigma, work->solver->enter, nb_active, TRUE, c);
  double one [2] = {1,0};
  QAtsigmaA = CHOLMOD(add)(work->data->Q, AtsigmaA, one, one, TRUE, FALSE, c);
  QAtsigmaA->stype = work->data->Q->stype;
  
  ldlchol(QAtsigmaA, work, c);

  CHOLMOD(free_sparse)(&AtsigmaA, c);
  CHOLMOD(free_sparse)(&QAtsigmaA, c);
}

void ldlupdate_entering_constraints(QPALMWorkspace *work, solver_common *c) {
  cholmod_sparse *Ae;
  Ae = CHOLMOD(submatrix)(work->solver->At_sqrt_sigma, NULL, -1, 
                      work->solver->enter, work->solver->nb_enter, TRUE, TRUE, c);
  //LD = ldlupdate(LD,Ae,'+');
  CHOLMOD(updown)(TRUE, Ae, work->solver->LD, c);
  CHOLMOD(free_sparse)(&Ae, c);
}

void ldldowndate_leaving_constraints(QPALMWorkspace *work, solver_common *c) {
  cholmod_sparse *Al;
  Al = CHOLMOD(submatrix)(work->solver->At_sqrt_sigma, NULL, -1, 
                      work->solver->leave, work->solver->nb_leave, TRUE, TRUE, c);
  //LD = ldlupdate(LD,Ae,'+');
  CHOLMOD(updown)(FALSE, Al, work->solver->LD, c);
  CHOLMOD(free_sparse)(&Al, c);
}

void ldlupdate_sigma_changed(QPALMWorkspace *work, solver_common *c) {
  cholmod_sparse *Ae;
  c_float *At_scalex = work->solver->At_scale->x;
  c_int *sigma_changed = work->solver->enter;

  size_t k, nb_sigma_changed = (size_t) work->nb_sigma_changed;
  for (k=0; k < nb_sigma_changed; k++) {
    At_scalex[sigma_changed[k]]= c_sqrt(1-1/(At_scalex[sigma_changed[k]]*At_scalex[sigma_changed[k]])); 
  }

  CHOLMOD(scale)(work->solver->At_scale, CHOLMOD_COL, work->solver->At_sqrt_sigma, c);
  Ae = CHOLMOD(submatrix)(work->solver->At_sqrt_sigma, NULL, -1, 
                      sigma_changed, work->nb_sigma_changed, TRUE, TRUE, c);
  for (k=0; k < work->data->m; k++) {
    At_scalex[k]=1.0/At_scalex[k]; 
  }
  CHOLMOD(scale)(work->solver->At_scale, CHOLMOD_COL, work->solver->At_sqrt_sigma, c);
  //LD = ldlupdate(LD,Ae,'+');
  CHOLMOD(updown)(TRUE, Ae, work->solver->LD, c);
  CHOLMOD(free_sparse)(&Ae, c);
}

void ldlsolveLD_neg_dphi(QPALMWorkspace *work, solver_common *c) {
  //set -dphi
  prea_vec_copy(work->dphi, work->neg_dphi, work->data->n);
  vec_self_mult_scalar(work->neg_dphi, -1, work->data->n);
  //d = ldlsolve(LD, -dphi)
  if (work->solver->d) {
      CHOLMOD(free_dense)(&work->solver->d, c);
  }
  work->solver->d = CHOLMOD(solve) (CHOLMOD_LDLt, work->solver->LD, work->solver->neg_dphi, c);
  work->d = work->solver->d->x;
}


void cholmod_set_settings(solver_common *c) {
  //Suitesparse memory allocation functions
  SuiteSparse_config.malloc_func = c_malloc;
  SuiteSparse_config.calloc_func = c_calloc;
  SuiteSparse_config.free_func = c_free;
  SuiteSparse_config.realloc_func = c_realloc;
  
  c->final_asis = FALSE ;
  c->final_super = FALSE ;
  c->final_ll = FALSE ;
  c->final_pack = TRUE ;
  c->final_monotonic = TRUE ;
  c->final_resymbol = TRUE ;
  c->quick_return_if_not_posdef = TRUE;
  c->nmethods = 1 ;
	c->method [0].ordering = CHOLMOD_NATURAL ;
	c->postorder = FALSE ;
  c->useGPU = FALSE ;
}

#endif /* USE_CHOLMOD */