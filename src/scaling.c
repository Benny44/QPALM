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
        // printf("Scaling\n");
        // for (size_t j = 0; j < work->data->n; j++) {
        //   printf("%f ", work->D_temp[j]);
        // }
        // printf("\n");
        // for (size_t j = 0; j < work->data->m; j++) {
        //   printf("%f ", work->E_temp[j]);
        // }
        // printf("\n");


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
        // printf("\n\n\n");
        // printf("Ax before scaling:\n");
        // c_float *Ax1 = work->data->A->x;
        // for (size_t j = 0; j < work->data->A->nzmax; j++) {
        // printf("%f ", Ax1[j]);
        // }
        // printf("\n");
        // printf("nzmax %d\n", (int)work->data->A->nzmax);
        // printf("nrow %d\n", (int)work->data->A->nrow);
        // printf("ncol %d\n", (int)work->data->A->ncol);
        // printf("itype %d\n", (int)work->data->A->itype);
        // printf("xtype %d\n", (int)work->data->A->xtype);
        // printf("dtype %d\n", (int)work->data->A->dtype);
        // printf("sorted %d\n", (int)work->data->A->sorted);
        // printf("packed %d\n", (int)work->data->A->packed);

        // printf("Ex: \n");
        // c_float *Ex = work->chol->E_temp->x;
        // for (size_t j = 0; j < work->chol->E_temp->nzmax; j++) {
        //   printf("%f ", Ex[j]);
        // }
        // printf("\n");
        // printf("E nrow: %d\n", (int)work->chol->E_temp->nrow);
        // printf("E ncol: %d\n", (int)work->chol->E_temp->ncol);
        // printf("E xtype: %d\n", (int)work->chol->E_temp->xtype);
        // printf("E dtype: %d\n", (int)work->chol->E_temp->dtype);
        // printf("E d: %d\n", (int)work->chol->E_temp->d);
        // // work->chol->E_temp->
        // printf("\n STATUS: %d \n",(int) (&work->chol->c)->status);
        
        CHOLMOD(scale)(work->chol->E_temp, CHOLMOD_ROW, work->data->A, &work->chol->c);

        // printf("\n STATUS: %d \n",(int) (&work->chol->c)->status);

        
        // printf("Ax after scaling with E:\n");
        // c_float *Ax = work->data->A->x;
        // for (size_t j = 0; j < work->data->A->nzmax; j++) {
        // printf("%f ", Ax[j]);
        // }
        // printf("\n");

        // printf("\n\n\n");

        // printf("Dx: \n");
        // c_float *Dx = work->chol->D_temp->x;
        // for (size_t j = 0; j < work->chol->D_temp->nzmax; j++) {
        //   printf("%f ", Dx[j]);
        // }
        // printf("\n");
        // printf("D nrow: %d\n", (int)work->chol->D_temp->nrow);
        // printf("D ncol: %d\n", (int)work->chol->D_temp->ncol);

        CHOLMOD(scale)(work->chol->D_temp, CHOLMOD_COL, work->data->A, &work->chol->c);

        // printf("Ax after scaling with D:\n");
        // c_float *Ax2 = work->data->A->x;
        // for (size_t j = 0; j < work->data->A->nzmax; j++) {
        // printf("%f ", Ax2[j]);
        // }
        // printf("\n");

        // printf("\n\n\n");


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