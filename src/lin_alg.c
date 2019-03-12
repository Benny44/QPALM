#include "lin_alg.h"


/* VECTOR FUNCTIONS ----------------------------------------------------------*/

c_float* vec_copy(c_float *a, c_int n) {
  c_float *b;
  c_int    i;

  b = c_malloc(n * sizeof(c_float));

  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }

  return b;
}

void prea_vec_copy(const c_float *a, c_float *b, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

void prea_int_vec_copy(const c_int *a, c_int *b, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

void vec_set_scalar(c_float *a, c_float sc, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    a[i] = sc;
  }
}

void vec_mult_scalar(c_float *a, c_float sc, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    a[i] *= sc;
  }
}

c_float vec_prod(const c_float *a, const c_float *b, c_int n) {
  c_float prod = 0.0;
  c_int   i; // Index

  for (i = 0; i <= n-4; i+=4) {
    prod += (a[i] * b[i] + a[i+1] * b[i+1] + a[i+2] * b[i+2] + a[i+3] * b[i+3]);
  }
  for (; i < n; i++) {
    prod += a[i] * b[i];
  }

  return prod;
}


void vec_ew_prod(const c_float *a, const c_float *b, c_float *c, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    c[i] = a[i] * b[i];
  }
}


void vec_ew_div(const c_float *a, const c_float *b, c_float *c, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    c[i] = a[i] / b[i];
  }
}


void vec_add_scaled(const c_float *a,
                    const c_float *b,
                    c_float       *c,
                    c_float        sc,
                    c_int          n) {
  c_int i;

  for (i = 0; i < n; i++) {
    c[i] =  a[i] + sc * b[i];
  }
}

c_float vec_norm_inf(const c_float *v, c_int l) {
  c_int   i;
  c_float abs_v_i;
  c_float max = 0.0;

  for (i = 0; i < l; i++) {
    abs_v_i = c_absval(v[i]);

    if (abs_v_i > max) max = abs_v_i;
  }
  return max;
}

void vec_ew_recipr(const c_float *a, c_float *b, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    b[i] = (c_float)1.0 / a[i];
  }
}

void vec_ew_max_vec(const c_float *a, const c_float *b, c_float *c, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    c[i] = c_max(a[i], b[i]);
  }
}

void vec_ew_min_vec(const c_float *a, const c_float *b, c_float *c, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    c[i] = c_min(a[i], b[i]);
  }
}

void vec_ew_mid_vec(const c_float *a, const c_float *bmin, const c_float *bmax, c_float *c, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    c[i] = c_max(bmin[i], c_min(a[i], bmax[i]));
  }
}

void vec_ew_sqrt(const c_float *a, c_float *b, c_int n) {
  c_int i;

  for (i = 0; i < n; i++) {
    b[i] = c_sqrt(a[i]);
  }
}



/* MATRIX FUNCTIONS ----------------------------------------------------------*/

void mat_premult_diag(csc *A, const c_float *d) {
  c_int j, i;

  for (j = 0; j < A->n; j++) {                // Cycle over columns
    for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row in the column
      A->x[i] *= d[A->i[i]];                  // Scale by corresponding element
                                              // of d for row i
    }
  }
}

void mat_postmult_diag(csc *A, const c_float *d) {
  c_int j, i;

  for (j = 0; j < A->n; j++) {                // Cycle over columns j
    for (i = A->p[j]; i < A->p[j + 1]; i++) { // Cycle every row i in column j
      A->x[i] *= d[j];                        // Scale by corresponding element
                                              // of d for column j
    }
  }
}


void mat_vec(const csc *A, const c_float *x, c_float *y) {
  c_int i, j;

  // y = 0
  for (i = 0; i < A->m; i++) {
    y[i] = 0;
  }

  // if A is empty
  if (A->p[A->n] == 0) {
    return;
  }

  // y =  A*x
  for (j = 0; j < A->n; j++) {
    for (i = A->p[j]; i < A->p[j + 1]; i++) {
      y[A->i[i]] += A->x[i] * x[j];
    }
  }
}

void mat_tpose_vec(const csc *A, const c_float *x, c_float *y) {
  c_int j, k;

  // if A is empty
  // if (A->p[A->n] == 0) {
  //   return;
  // }
  // y =  A*x
  for (j = 0; j < A->n; j++) {
    y[j] = 0;
    for (k = A->p[j]; k <= A->p[j + 1]-4; k+=4) {
      y[j] += (A->x[k] * x[A->i[k]] 
              + A->x[k+1] * x[A->i[k+1]] 
              + A->x[k+2] * x[A->i[k+2]] 
              + A->x[k+3] * x[A->i[k+3]]);
    }
    for(; k < A->p[j + 1]; k++) {
        y[j] += A->x[k] * x[A->i[k]];
    }
    // for (k = A->p[j]; k < A->p[j + 1]; k++) {
    //   y[j] += A->x[k] * x[A->i[k]];
    // }
  }
}

void mat_tpose_vec_skip_diag_add(const csc *A, const c_float *x, c_float *y) {
  c_int i, j, k;

  // if A is empty
  if (A->p[A->n] == 0) {
    return;
  }
  // y +=  A*x
    for (j = 0; j < A->n; j++) {
      for (k = A->p[j]; k < A->p[j + 1]; k++) {
        i     = A->i[k];
        y[j] += i == j ? 0 : A->x[k] * x[i];
      }
    }
}

void mat_vec_triu(const csc *Q, const c_float *x, c_float *y) {
  // Q * x (upper triangular part)
  mat_vec(Q, x, y);

  // Q' * x (lower triangular part with no diagonal)
  mat_tpose_vec_skip_diag_add(Q, x, y);
}


void mat_inf_norm_cols(const csc *M, c_float *E) {
  c_int j, ptr;

  // Initialize zero max elements
  for (j = 0; j < M->n; j++) {
    E[j] = 0.;
  }

  // Compute maximum across columns
  for (j = 0; j < M->n; j++) {
    for (ptr = M->p[j]; ptr < M->p[j + 1]; ptr++) {
      E[j] = c_max(c_absval(M->x[ptr]), E[j]);
    }
  }
}

void mat_inf_norm_rows(const csc *M, c_float *E) {
  c_int i, j, ptr;

  // Initialize zero max elements
  for (j = 0; j < M->m; j++) {
    E[j] = 0.;
  }

  // Compute maximum across rows
  for (j = 0; j < M->n; j++) {
    for (ptr = M->p[j]; ptr < M->p[j + 1]; ptr++) {
      i    = M->i[ptr];
      E[i] = c_max(c_absval(M->x[ptr]), E[i]);
    }
  }
}

