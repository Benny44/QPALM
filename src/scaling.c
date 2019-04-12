/**
 * @file scaling.c
 * @author Ben Hermans
 * @brief Problem data scaling during setup.
 * @details This file includes the routine that is called during setup to scale the problem data, 
 * and initial guesses if the problem is warm-started. Scaling the problem is useful to prevent 
 * large changes in the active set and to guard against ill-conditioning in the objective function. 
 * @note The function in this file makes use of the cholmod scale routines. 
 */
#include "scaling.h"
#include "cholmod.h"
#include <stdio.h>

// Set values lower than threshold MIN_SCALING to 1 and larger than MAX_SCALING to MAX_SCALING
void limit_scaling(c_float *D, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    D[i] = D[i] < MIN_SCALING ? 1.0 : D[i];
    D[i] = D[i] > MAX_SCALING ? MAX_SCALING : D[i];
  }
}

void scale_data(QPALMWorkspace *work) {

    vec_set_scalar(work->scaling->D, 1, work->data->n);
    vec_set_scalar(work->scaling->Dinv, 1, work->data->n);
    vec_set_scalar(work->scaling->E, 1, work->data->m);
    vec_set_scalar(work->scaling->Einv, 1, work->data->m);
    
    c_int i;
    //Ruiz on constraint matrix A
    for (i = 0; i < work->settings->scaling; i++) {

        // Set D_temp = vecnorm(A,inf,1) (cols) and E_temp = vecnorm(A,inf,2) (rows)
        mat_inf_norm_cols(work->data->A, work->D_temp);
        mat_inf_norm_rows(work->data->A, work->E_temp);

        // Set to 1 values with 0 norms (avoid crazy scaling)
        limit_scaling(work->D_temp, work->data->n);
        limit_scaling(work->E_temp, work->data->m);

        // Take square root of norms
        vec_ew_sqrt(work->D_temp, work->D_temp, work->data->n);
        vec_ew_sqrt(work->E_temp, work->E_temp, work->data->m);

        // 1./D and 1./E
        vec_ew_recipr(work->D_temp, work->D_temp, work->data->n);
        vec_ew_recipr(work->E_temp, work->E_temp, work->data->m);

        // Equilibrate matrix A
        // A <- EAD
        CHOLMOD(scale)(work->chol->E_temp, CHOLMOD_ROW, work->data->A, &work->chol->c);
        CHOLMOD(scale)(work->chol->D_temp, CHOLMOD_COL, work->data->A, &work->chol->c);

        // Update equilibration matrices D and E
        vec_ew_prod(work->scaling->D, work->D_temp, work->scaling->D, work->data->n);
        vec_ew_prod(work->scaling->E, work->E_temp, work->scaling->E, work->data->m);

    }

    // Equilibrate matrix Q and vector q
    // Q <- DQD, q <- Dq
    prea_vec_copy(work->scaling->D, work->D_temp, work->data->n);
    CHOLMOD(scale)(work->chol->D_temp, CHOLMOD_SYM, work->data->Q, &work->chol->c);
    vec_ew_prod(work->scaling->D, work->data->q, work->data->q, work->data->n);

    // Cost scaling
    vec_add_scaled(work->Qx, work->data->q, work->temp_n, 1, work->data->n);
    work->scaling->c = 1/c_max(1.0, vec_norm_inf(work->temp_n, work->data->n));
    vec_mult_scalar(work->data->q, work->scaling->c, work->data->n);
    cholmod_dense *c = CHOLMOD(ones)(1,1,CHOLMOD_REAL, &work->chol->c);
    c_float *cx = c->x;
    cx[0] = work->scaling->c;
    CHOLMOD(scale)(c, CHOLMOD_SCALAR, work->data->Q, &work->chol->c);
    CHOLMOD(free_dense)(&c, &work->chol->c);

    // Store cinv, Dinv, Einv
    vec_ew_recipr(work->scaling->D, work->scaling->Dinv, work->data->n);
    vec_ew_recipr(work->scaling->E, work->scaling->Einv, work->data->m);
    work->scaling->cinv = 1/work->scaling->c;

    // Scale problem vectors l, u
    vec_ew_prod(work->scaling->E, work->data->bmin, work->data->bmin, work->data->m);
    vec_ew_prod(work->scaling->E, work->data->bmax, work->data->bmax, work->data->m);

    // Scale initial vectors x, xprev, x0, Qx, Ax, and y if they are warm-started
    if (work->settings->warm_start) {
      vec_ew_prod(work->x, work->scaling->Dinv, work->x, work->data->n);
      prea_vec_copy(work->x, work->x0, work->data->n);
      prea_vec_copy(work->x, work->x_prev, work->data->n);

      vec_ew_prod(work->Qx, work->scaling->D, work->Qx, work->data->n);
      vec_mult_scalar(work->Qx, work->scaling->c, work->data->n);

      vec_ew_prod(work->Ax,work->scaling->E, work->Ax, work->data->m);

      vec_ew_prod(work->y, work->scaling->Einv, work->y, work->data->m);
      vec_mult_scalar(work->y, work->scaling->c, work->data->m);
    }  

}