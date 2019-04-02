#include "cholmod_interface.h"
#include <stdio.h>

void mat_vec(cholmod_sparse *A, cholmod_dense *x, cholmod_dense *y, cholmod_common *c) {
    double one [2] = {1,0};
    cholmod_sdmult(A, 0, one, one, x, y, c);
}

void mat_tpose_vec(cholmod_sparse *A, cholmod_dense *x, cholmod_dense *y, cholmod_common *c) {
    double one [2] = {1,0};
    cholmod_sdmult(A, 1, one, one, x, y, c);
}

void mat_inf_norm_cols(cholmod_sparse *M, c_float *E) {
  int j, i;
  c_float *Mx = M->x;
  int *Mp = M->p;

  // Initialize zero max elements
  for (j = 0; j < M->ncol; j++) {
    E[j] = 0.;
  }

  // Compute maximum across columns
  for (j = 0; j < M->ncol; j++) {
    for (i = Mp[j]; i < Mp[j + 1]; i++) {
      E[j] = c_max(c_absval(Mx[i]), E[j]);
    }
  }
}

void mat_inf_norm_rows(cholmod_sparse *M, c_float *E) {
  int i, j, k;
  int *Mp = M->p;
  int *Mi = M->i;
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