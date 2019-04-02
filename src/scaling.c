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

        // // Set to 1 values with 0 norms (avoid crazy scaling)
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
        // mat_premult_diag(work->data->A, work->E_temp);
        // mat_postmult_diag(work->data->A, work->D_temp);
        cholmod_scale(work->chol->E_temp, CHOLMOD_ROW, work->data->A, &work->chol->c);
        cholmod_scale(work->chol->D_temp, CHOLMOD_COL, work->data->A, &work->chol->c);

        // Update equilibration matrices D and E
        vec_ew_prod(work->scaling->D, work->D_temp, work->scaling->D, work->data->n);
        vec_ew_prod(work->scaling->E, work->E_temp, work->scaling->E, work->data->m);

    }

    // Equilibrate matrix Q and vector q
    // Q <- DPD, q <- Dq
    prea_vec_copy(work->scaling->D, work->D_temp, work->data->n);
    cholmod_scale(work->chol->D_temp, CHOLMOD_SYM, work->data->Q, &work->chol->c);
    // mat_premult_diag(work->data->Q, work->scaling->D);
    // mat_postmult_diag(work->data->Q, work->scaling->D);
    vec_ew_prod(work->scaling->D, work->data->q, work->data->q, work->data->n);

    // Store cinv, Dinv, Einv
    vec_ew_recipr(work->scaling->D, work->scaling->Dinv, work->data->n);
    vec_ew_recipr(work->scaling->E, work->scaling->Einv, work->data->m);


    // Scale problem vectors l, u
    vec_ew_prod(work->scaling->E, work->data->bmin, work->data->bmin, work->data->m);
    vec_ew_prod(work->scaling->E, work->data->bmax, work->data->bmax, work->data->m);

}