/**
 * @file scaling.c
 * @author Ben Hermans
 * @brief Problem data scaling during setup.
 * @details This file includes the routine that is called during setup to scale the problem data, 
 * and initial guesses if the problem is warm-started. Scaling the problem is useful to prevent 
 * large changes in the active set and to guard against ill-conditioning in the objective function. 
 * @note The function in this file makes use of the cholmod scale routines. 
 */


# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include <stdio.h>
#include "scaling.h"
#include "lin_alg.h"
#ifdef USE_LADEL
#include "ladel.h"
#elif defined USE_CHOLMOD
#include "cholmod.h"
#endif

// Set values lower than threshold MIN_SCALING to 1
void limit_scaling(c_float *D, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    D[i] = D[i] < MIN_SCALING ? 1.0 : D[i];
  }
}

void scale_data(QPALMWorkspace *work) {
    solver_common c;
    #ifdef USE_LADEL
    #elif defined USE_CHOLMOD
    CHOLMOD(start)(&c);
    cholmod_set_settings(&c);
    #endif

    size_t n = work->data->n;
    size_t m = work->data->m;
    vec_set_scalar(work->scaling->D, 1, n);
    vec_set_scalar(work->scaling->E, 1, m);
    
    c_int i;
    //Ruiz on constraint matrix A
    for (i = 0; i < work->settings->scaling; i++) {

        // Set D_temp = vecnorm(A,inf,1) (cols) and E_temp = vecnorm(A,inf,2) (rows)
        mat_inf_norm_cols(work->data->A, work->D_temp);
        mat_inf_norm_rows(work->data->A, work->E_temp);

        // Set to 1 values with 0 norms (avoid crazy scaling)
        limit_scaling(work->D_temp, n);
        limit_scaling(work->E_temp, m);

        // Take square root of norms
        vec_ew_sqrt(work->D_temp, work->D_temp, n);
        vec_ew_sqrt(work->E_temp, work->E_temp, m);

        // 1./D and 1./E
        vec_ew_recipr(work->D_temp, work->D_temp, n);
        vec_ew_recipr(work->E_temp, work->E_temp, m);

        // Equilibrate matrix A
        // A <- EAD
        
        #ifdef USE_LADEL
        ladel_scale_rows(work->data->A, work->solver->E_temp);
        ladel_scale_columns(work->data->A, work->solver->D_temp);
        #elif defined USE_CHOLMOD
        CHOLMOD(scale)(work->solver->E_temp, CHOLMOD_ROW, work->data->A, &c);
        CHOLMOD(scale)(work->solver->D_temp, CHOLMOD_COL, work->data->A, &c);
        #endif
        // Update equilibration matrices D and E
        vec_ew_prod(work->scaling->D, work->D_temp, work->scaling->D, n);
        vec_ew_prod(work->scaling->E, work->E_temp, work->scaling->E, m);
    }

    // Equilibrate matrix Q and vector q
    // Q <- cDQD, q <- cDq
    vec_ew_prod(work->scaling->D, work->data->q, work->data->q, n);
    vec_ew_prod(work->scaling->D, work->Qx, work->Qx, n);
    prea_vec_copy(work->scaling->D, work->D_temp, n);
    vec_add_scaled(work->Qx, work->data->q, work->dphi, 1, n);
    work->scaling->c = 1/c_max(1.0, vec_norm_inf(work->dphi, n));
    vec_self_mult_scalar(work->data->q, work->scaling->c, n);
    
    #ifdef USE_LADEL
    ladel_scale_columns(work->data->Q, work->solver->D_temp);
    ladel_scale_rows(work->data->Q, work->solver->D_temp);
    ladel_scale_scalar(work->data->Q, work->scaling->c);
    #elif defined USE_CHOLMOD
    CHOLMOD(scale)(work->solver->D_temp, CHOLMOD_SYM, work->data->Q, &c);
    cholmod_dense *scalar = CHOLMOD(ones)(1,1,CHOLMOD_REAL, &c);
    c_float *scalarx = scalar->x;
    scalarx[0] = work->scaling->c;
    CHOLMOD(scale)(scalar, CHOLMOD_SCALAR, work->data->Q, &c);
    CHOLMOD(free_dense)(&scalar, &c);
    CHOLMOD(finish)(&c);
    #endif

    // Store cinv, Dinv, Einv
    vec_ew_recipr(work->scaling->D, work->scaling->Dinv, n);
    vec_ew_recipr(work->scaling->E, work->scaling->Einv, m);
    work->scaling->cinv = (c_float) 1.0/work->scaling->c;

    // Scale problem vectors l, u
    vec_ew_prod(work->scaling->E, work->data->bmin, work->data->bmin, m);
    vec_ew_prod(work->scaling->E, work->data->bmax, work->data->bmax, m);    
}


# ifdef __cplusplus
}
# endif // ifdef __cplusplus