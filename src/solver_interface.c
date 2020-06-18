/**
 * @file solver_interface.c
 * @author Ben Hermans
 * @brief Interface and wrapper to matrix/factorization (ladel/cholmod) functions
 * @details This file includes all calls to cholmod/ladel functions apart from scaling in scaling.c and memory
 * allocation/deallocation in the main functions in qpalm.c. It includes all matrix operations, such as
 * matrix vector products, row- and columnwise norms, cholesky factorizations, factorization updates and
 * solving the linear system. Finally, all the settings relevant to cholmod (and suitesparse) are included
 * in this file as well.
 */

#include "lin_alg.h"
#include "solver_interface.h"
#include <stdio.h>

#ifdef USE_LADEL
#include "ladel.h"
#endif

void qpalm_set_factorization_method(QPALMWorkspace *work)
{
  #ifdef USE_LADEL
  if (work->settings->factorization_method == FACTORIZE_KKT_OR_SCHUR)
  {
    /* TODO: determine criterion to set the factorization method depending on Q and A */
    work->solver->factorization_method = FACTORIZE_KKT;
  } else
  {
    work->solver->factorization_method = work->settings->factorization_method;
  }

  #elif defined USE_CHOLMOD
  work->solver->factorization_method = FACTORIZE_SCHUR;
  #endif
}


#ifdef USE_LADEL

void mat_vec(solver_sparse *A, solver_dense *x, solver_dense *y, solver_common *c) 
{
    ladel_int n = A->ncol;
    if (x!=y) {
      if (A->symmetry == UNSYMMETRIC)
        ladel_matvec(A, x, y, TRUE);
      else
        ladel_symmetric_matvec(A, x, y, TRUE);
    } else {
      ladel_double* x2 = ladel_malloc(n, sizeof(c_float));
      ladel_double_vector_copy(x, n, x2);
      if (A->symmetry == UNSYMMETRIC)
        ladel_matvec(A, x2, y, TRUE);
      else
        ladel_symmetric_matvec(A, x2, y, TRUE);
      ladel_free(x2);
    }
}

void mat_tpose_vec(solver_sparse *A, solver_dense *x, solver_dense *y, solver_common *c)
{
    ladel_int m = A->nrow;
    if (x!=y) {
      if (A->symmetry == UNSYMMETRIC)
        ladel_tpose_matvec(A, x, y, TRUE);
      else
        ladel_symmetric_matvec(A, x, y, TRUE);  
    } else {
      ladel_double* x2 = ladel_malloc(m, sizeof(c_float));
      ladel_double_vector_copy(x, m, x2);
      if (A->symmetry == UNSYMMETRIC)
        ladel_tpose_matvec(A, x2, y, TRUE);
      else
        ladel_symmetric_matvec(A, x2, y, TRUE); 
      ladel_free(x2);
    }
}

// TODO: incorporate the factorizations here also, and write a version of FACTORIZE_KKT using cholmod 
void qpalm_form_kkt(QPALMWorkspace *work)
{
    solver_sparse *Q = work->data->Q, *A = work->data->A, *kkt = work->solver->kkt, *kkt_full = work->solver->kkt_full, *At = work->solver->At;
    ladel_int col, index, index_kkt, n = work->data->n, m = work->data->m, Qnz = Q->nzmax;
    c_float *sigma_inv = work->sigma_inv, *first_elem_A = work->solver->first_elem_A;
    c_int *first_row_A = work->solver->first_row_A;
    /* copy Q */
    for (col = 0; col < n; col++)
    {
        kkt->p[col] = kkt_full->p[col] = Q->p[col];
        kkt->nz[col] = Q->p[col+1] - Q->p[col];
    }
    kkt->p[col] = kkt_full->p[col] = Q->p[col];
    for (index = 0; index < Qnz; index++)
    {
        kkt->i[index] = kkt_full->i[index] = Q->i[index];
        kkt->x[index] = kkt_full->x[index] = Q->x[index];
    }

    /* copy [At; -\Sigma^{-1}] */
    index_kkt = Qnz;
    for (; col < m+n; col++)
    {
        kkt_full->i[index_kkt] = first_row_A[col-n] = At->i[At->p[col-n]];
        kkt_full->x[index_kkt] = first_elem_A[col-n] = At->x[At->p[col-n]];

        if (work->solver->active_constraints[col-n])
        {
            kkt->nz[col] = At->p[col-n+1] - At->p[col-n] + 1;
            kkt->i[index_kkt] = At->i[At->p[col-n]];
            kkt->x[index_kkt] = At->x[At->p[col-n]];
        } 
        else 
        {
            kkt->nz[col] = 1;
            kkt->i[index_kkt] = col;
            kkt->x[index_kkt] = 1;
        }

        if (At->p[col-n+1]-At->p[col-n] != 0) index_kkt++;

        for (index = At->p[col-n]+1; index < At->p[col-n+1]; index++)
        {
            kkt->i[index_kkt] = kkt_full->i[index_kkt] = At->i[index];
            kkt->x[index_kkt] = kkt_full->x[index_kkt] = At->x[index];
            index_kkt++;
        }

        kkt->i[index_kkt] = kkt_full->i[index_kkt] = col;
        kkt->x[index_kkt] = kkt_full->x[index_kkt] = -sigma_inv[col-n];
        if (At->p[col-n+1]-At->p[col-n] == 0) kkt->x[index_kkt] = 1;
        index_kkt++;

        kkt->p[col+1] = kkt_full->p[col+1] = Qnz + At->p[col+1-n] + 1 + col - n;
    }
}


void qpalm_reform_kkt(QPALMWorkspace *work)
{
    solver_sparse *kkt = work->solver->kkt, *At = work->solver->At;
    ladel_int col, index, n = work->data->n, m = work->data->m;
    ladel_int *first_row_A = work->solver->first_row_A;
    ladel_double *sigma_inv = work->sigma_inv, *first_elem_A = work->solver->first_elem_A;

    for (col = n; col < n+m; col++)
    {
        if (work->solver->active_constraints[col-n])
        {
            kkt->nz[col] = At->p[col-n+1] - At->p[col-n] + 1;
            kkt->i[kkt->p[col]] = first_row_A[col-n];
            kkt->x[kkt->p[col]] = first_elem_A[col-n];
            kkt->x[kkt->p[col+1]-1] = -sigma_inv[col-n]; /* row index should be correct already */
            kkt->i[kkt->p[col+1]-1] = col;
        } else
        {
            kkt->nz[col] = 1;
            kkt->i[kkt->p[col]] = col;
            kkt->x[kkt->p[col]] = 1; 
        }
    }
}

void kkt_update_entering_constraints(QPALMWorkspace *work, solver_common *c)
{
    solver_sparse *kkt = work->solver->kkt, *At = work->solver->At;
    ladel_int col, index, n = work->data->n, m = work->data->m;
    ladel_int *first_row_A = work->solver->first_row_A;
    ladel_double *sigma_inv = work->sigma_inv, *first_elem_A = work->solver->first_elem_A;

    for (index = 0; index < work->solver->nb_enter; index++)
    {
        col = work->solver->enter[index] + n;
        kkt->nz[col] = At->p[col-n+1] - At->p[col-n] + 1;
        kkt->i[kkt->p[col]] = first_row_A[col-n];
        kkt->x[kkt->p[col]] = first_elem_A[col-n];
        kkt->x[kkt->p[col+1]-1] = -sigma_inv[col-n]; /* row index should be correct already */
        ladel_row_add(work->solver->LD, work->solver->sym, col, kkt, col, -sigma_inv[col-n], c);
    }
}

void kkt_update_leaving_constraints(QPALMWorkspace *work, solver_common *c)
{
    ladel_int col, index, n = work->data->n, m = work->data->m;

    for (index = 0; index < work->solver->nb_leave; index++)
    {
        col = work->solver->leave[index] + n;
        ladel_row_del(work->solver->LD, work->solver->sym, col, c);
        /* no need to update the kkt system here */
    }
}

void kkt_solve(QPALMWorkspace *work, solver_common *c)
{
    c_int n = work->data->n, m = work->data->m;
    prea_vec_copy(work->dphi, work->solver->rhs_kkt, n);
    vec_self_mult_scalar(work->solver->rhs_kkt, -1, n);
    vec_set_scalar(work->solver->rhs_kkt + n, 0, m);

    ladel_dense_solve(work->solver->LD, work->solver->rhs_kkt, work->solver->sol_kkt, c);
    prea_vec_copy(work->solver->sol_kkt, work->d, n);
}

#elif defined USE_CHOLMOD
#include "cholmod.h"

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

#endif /* USE_CHOLMOD */


void ldlchol(solver_sparse *M, QPALMWorkspace *work, solver_common *c) {
  #ifdef USE_LADEL
  ladel_diag d;
  d.diag_elem = 1.0/work->gamma;
  if (work->settings->proximal) d.diag_size = work->data->n;
  else d.diag_size = 0;

  if (work->solver->first_factorization)
  {
    work->solver->LD = ladel_factor_free(work->solver->LD);
    /* Compute the pattern of Q+A^T*A to allocate L */
    solver_sparse *AtA, *QAtA;
    AtA = ladel_mat_mat_transpose_pattern(work->solver->At_sqrt_sigma, work->data->A, c);
    QAtA = ladel_add_matrices_pattern(work->data->Q, AtA, c);
    QAtA->symmetry = UPPER;

    /* TODO: consider SCHUR method also with ordering */
    ladel_factorize_advanced_with_diag(M, d, work->solver->sym, NO_ORDERING, &work->solver->LD, QAtA, c);
    
    ladel_sparse_free(AtA);
    ladel_sparse_free(QAtA);
    work->solver->first_factorization = FALSE;
  }
  else
  {
    ladel_factorize_with_prior_basis_with_diag(M, d, work->solver->sym, work->solver->LD, c);
  }
  #elif defined USE_CHOLMOD
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
  #endif /* USE_CHOLMOD */
}

void ldlcholQAtsigmaA(QPALMWorkspace *work, solver_common *c) {
  solver_sparse *AtsigmaA;
  solver_sparse *QAtsigmaA;
  size_t nb_active = 0;
  for (size_t i = 0; i < work->data->m; i++) {
      if (work->solver->active_constraints[i]){
          work->solver->enter[nb_active] = (c_int)i;
          nb_active++;
      }      
  }
  #ifdef USE_LADEL
  solver_sparse *At_sqrt_sigma = ladel_column_submatrix(work->solver->At_sqrt_sigma, work->solver->enter, nb_active);
  solver_sparse *A_sqrt_sigma = ladel_transpose(At_sqrt_sigma, TRUE, c);
  AtsigmaA = ladel_mat_mat_transpose(At_sqrt_sigma, A_sqrt_sigma, c);
  QAtsigmaA = ladel_add_matrices(1.0, work->data->Q, 1.0, AtsigmaA, c);
  QAtsigmaA->symmetry = UPPER;
  #elif defined USE_CHOLMOD
  AtsigmaA = CHOLMOD(aat)(work->solver->At_sqrt_sigma, work->solver->enter, nb_active, TRUE, c);
  double one [2] = {1,0};
  QAtsigmaA = CHOLMOD(add)(work->data->Q, AtsigmaA, one, one, TRUE, FALSE, c);
  QAtsigmaA->stype = work->data->Q->stype;
  #endif
  ldlchol(QAtsigmaA, work, c);

  #ifdef USE_LADEL
  AtsigmaA = ladel_sparse_free(AtsigmaA);
  QAtsigmaA = ladel_sparse_free(QAtsigmaA);
  At_sqrt_sigma = ladel_sparse_free(At_sqrt_sigma);
  A_sqrt_sigma = ladel_sparse_free(A_sqrt_sigma);
  #elif defined USE_CHOLMOD
  CHOLMOD(free_sparse)(&AtsigmaA, c);
  CHOLMOD(free_sparse)(&QAtsigmaA, c);
  #endif
}

void ldlupdate_entering_constraints(QPALMWorkspace *work, solver_common *c) {
  #ifdef USE_LADEL
  ladel_int index;
  for (index = 0; index < work->solver->nb_enter; index++)
  {
    ladel_rank1_update(work->solver->LD, work->solver->sym, work->solver->At_sqrt_sigma, 
                        work->solver->enter[index], 1.0, UPDATE, c);
  }
  #elif defined USE_CHOLMOD
  solver_sparse *Ae;
  Ae = CHOLMOD(submatrix)(work->solver->At_sqrt_sigma, NULL, -1, 
                      work->solver->enter, work->solver->nb_enter, TRUE, TRUE, c);
  //LD = ldlupdate(LD,Ae,'+');
  CHOLMOD(updown)(TRUE, Ae, work->solver->LD, c);
  CHOLMOD(free_sparse)(&Ae, c);
  #endif
}

void ldldowndate_leaving_constraints(QPALMWorkspace *work, solver_common *c) {
  #ifdef USE_LADEL
  ladel_int index;
  for (index = 0; index < work->solver->nb_leave; index++)
  {
    ladel_rank1_update(work->solver->LD, work->solver->sym, work->solver->At_sqrt_sigma, 
                        work->solver->leave[index], 1.0, DOWNDATE, c);
  }
  #elif defined USE_CHOLMOD
  cholmod_sparse *Al;
  Al = CHOLMOD(submatrix)(work->solver->At_sqrt_sigma, NULL, -1, 
                      work->solver->leave, work->solver->nb_leave, TRUE, TRUE, c);
  //LD = ldlupdate(LD,Ae,'+');
  CHOLMOD(updown)(FALSE, Al, work->solver->LD, c);
  CHOLMOD(free_sparse)(&Al, c);
  #endif
}

void ldlupdate_sigma_changed(QPALMWorkspace *work, solver_common *c) {
  c_int *sigma_changed = work->solver->enter;
  size_t k, nb_sigma_changed = (size_t) work->nb_sigma_changed;
  
  #ifdef USE_LADEL
  c_float *At_scalex = work->solver->At_scale;
  #elif defined USE_CHOLMOD
  solver_sparse *Ae;
  c_float *At_scalex = work->solver->At_scale->x;
  #endif
  
  for (k = 0; k < nb_sigma_changed; k++) {
    At_scalex[sigma_changed[k]]= c_sqrt(1-1/(At_scalex[sigma_changed[k]]*At_scalex[sigma_changed[k]])); 
  }

  #ifdef USE_LADEL
  for (k = 0; k < nb_sigma_changed; k++)
  {
    ladel_rank1_update(work->solver->LD, work->solver->sym, work->solver->At_sqrt_sigma, 
                        sigma_changed[k], At_scalex[sigma_changed[k]], UPDATE, c);
  }
  #elif defined USE_CHOLMOD
  CHOLMOD(scale)(work->solver->At_scale, CHOLMOD_COL, work->solver->At_sqrt_sigma, c);
  Ae = CHOLMOD(submatrix)(work->solver->At_sqrt_sigma, NULL, -1, 
                      sigma_changed, work->nb_sigma_changed, TRUE, TRUE, c);
  //LD = ldlupdate(LD,Ae,'+');
  CHOLMOD(updown)(TRUE, Ae, work->solver->LD, c);
  CHOLMOD(free_sparse)(&Ae, c);

  for (k = 0; k < work->data->m; k++) {
    At_scalex[sigma_changed[k]]= 1.0/At_scalex[sigma_changed[k]]; 
  }
  CHOLMOD(scale)(work->solver->At_scale, CHOLMOD_COL, work->solver->At_sqrt_sigma, c);
  #endif  
}

void ldlsolveLD_neg_dphi(QPALMWorkspace *work, solver_common *c) {
  //set -dphi
  prea_vec_copy(work->dphi, work->neg_dphi, work->data->n);
  vec_self_mult_scalar(work->neg_dphi, -1, work->data->n);
  #ifdef USE_LADEL
  ladel_dense_solve(work->solver->LD, work->neg_dphi, work->d, c);
  #elif defined USE_CHOLMOD
  //d = ldlsolve(LD, -dphi)
  if (work->solver->d) {
      CHOLMOD(free_dense)(&work->solver->d, c);
  }
  work->solver->d = CHOLMOD(solve) (CHOLMOD_LDLt, work->solver->LD, work->solver->neg_dphi, c);
  work->d = work->solver->d->x;
  #endif
}

#ifdef USE_LADEL
#elif defined USE_CHOLMOD
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