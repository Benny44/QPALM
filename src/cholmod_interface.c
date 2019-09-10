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

#include "cholmod_interface.h"
#include "lin_alg.h"
#include <stdio.h>

void mat_vec(cholmod_sparse *A, cholmod_dense *x, cholmod_dense *y, cholmod_common *c) {
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

void mat_tpose_vec(cholmod_sparse *A, cholmod_dense *x, cholmod_dense *y, cholmod_common *c) {
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

void mat_inf_norm_cols(cholmod_sparse *M, c_float *E) {
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

void mat_inf_norm_rows(cholmod_sparse *M, c_float *E) {
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

void ldlchol(cholmod_sparse *M, QPALMWorkspace *work) {
  if (work->chol->LD) {
      CHOLMOD(free_factor)(&work->chol->LD, &work->chol->c);
  }
  work->chol->LD = CHOLMOD(analyze) (M, &work->chol->c) ;
  if (work->settings->proximal) {
    double beta [2] = {1.0/work->gamma,0};
    CHOLMOD(factorize_p) (M, beta, NULL, 0, work->chol->LD, &work->chol->c);
  } else {
    CHOLMOD(factorize) (M, work->chol->LD, &work->chol->c);
  }
  if ((&work->chol->c)->status != CHOLMOD_OK) {
      (&work->chol->c)->supernodal = CHOLMOD_SIMPLICIAL;
      if (work->settings->proximal) {
        double beta [2] = {1.0/work->gamma,0};
        CHOLMOD(factorize_p) (M, beta, NULL, 0, work->chol->LD, &work->chol->c);
      } else {
        CHOLMOD(factorize) (M, work->chol->LD, &work->chol->c);
      }
  }

}

void ldlcholQ(QPALMWorkspace *work) {
  ldlchol(work->data->Q, work);
}

void ldlcholQAtsigmaA(QPALMWorkspace *work) {
  cholmod_sparse *AtsigmaA;
  cholmod_sparse *QAtsigmaA;
  size_t nb_active = 0;
  for (size_t i = 0; i < work->data->m; i++) {
      if (work->chol->active_constraints[i]){
          work->chol->enter[nb_active] = i;
          nb_active++;
      }      
  }
  AtsigmaA = CHOLMOD(aat)(work->chol->At_sqrt_sigma, work->chol->enter, nb_active, TRUE, &work->chol->c);
  double one [2] = {1,0};
  QAtsigmaA = CHOLMOD(add)(work->data->Q, AtsigmaA, one, one, TRUE, FALSE, &work->chol->c);
  QAtsigmaA->stype = work->data->Q->stype;
  
  ldlchol(QAtsigmaA, work);

  CHOLMOD(free_sparse)(&AtsigmaA, &work->chol->c);
  CHOLMOD(free_sparse)(&QAtsigmaA, &work->chol->c);
}

void ldlupdate_entering_constraints(QPALMWorkspace *work) {
  cholmod_sparse *Ae;
  Ae = CHOLMOD(submatrix)(work->chol->At_sqrt_sigma, NULL, -1, 
                      work->chol->enter, work->chol->nb_enter, TRUE, TRUE, &work->chol->c);
  //LD = ldlupdate(LD,Ae,'+');
  CHOLMOD(updown)(TRUE, Ae, work->chol->LD, &work->chol->c);
  CHOLMOD(free_sparse)(&Ae, &work->chol->c);
}

void ldldowndate_leaving_constraints(QPALMWorkspace *work) {
  cholmod_sparse *Al;
  Al = CHOLMOD(submatrix)(work->chol->At_sqrt_sigma, NULL, -1, 
                      work->chol->leave, work->chol->nb_leave, TRUE, TRUE, &work->chol->c);
  //LD = ldlupdate(LD,Ae,'+');
  CHOLMOD(updown)(FALSE, Al, work->chol->LD, &work->chol->c);
  CHOLMOD(free_sparse)(&Al, &work->chol->c);
}

void ldlupdate_sigma_changed(QPALMWorkspace *work) {
  cholmod_sparse *Ae;
  c_float *At_scalex = work->chol->At_scale->x;
  c_int *sigma_changed = work->chol->enter;

  size_t k;
  for (k=0; k < work->nb_sigma_changed; k++) {
    At_scalex[sigma_changed[k]]= c_sqrt(1-1/(At_scalex[sigma_changed[k]]*At_scalex[sigma_changed[k]])); 
  }

  CHOLMOD(scale)(work->chol->At_scale, CHOLMOD_COL, work->chol->At_sqrt_sigma, &work->chol->c);
  Ae = CHOLMOD(submatrix)(work->chol->At_sqrt_sigma, NULL, -1, 
                      sigma_changed, work->nb_sigma_changed, TRUE, TRUE, &work->chol->c);
  for (k=0; k < work->data->m; k++) {
    At_scalex[k]=1.0/At_scalex[k]; 
  }
  CHOLMOD(scale)(work->chol->At_scale, CHOLMOD_COL, work->chol->At_sqrt_sigma, &work->chol->c);
  //LD = ldlupdate(LD,Ae,'+');
  CHOLMOD(updown)(TRUE, Ae, work->chol->LD, &work->chol->c);
  CHOLMOD(free_sparse)(&Ae, &work->chol->c);
}

void ldlsolveLD_neg_dphi(QPALMWorkspace *work) {
  //set -dphi
  prea_vec_copy(work->dphi, work->neg_dphi, work->data->n);
  vec_mult_scalar(work->neg_dphi, -1, work->data->n);
  //d = ldlsolve(LD, -dphi)
  if (work->chol->d) {
      CHOLMOD(free_dense)(&work->chol->d, &work->chol->c);
  }
  work->chol->d = CHOLMOD(solve) (CHOLMOD_LDLt, work->chol->LD, work->chol->neg_dphi, &work->chol->c);
  work->d = work->chol->d->x;
}


void cholmod_set_settings(cholmod_common *c) {
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