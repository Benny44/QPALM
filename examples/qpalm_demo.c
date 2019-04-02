#include "qpalm.h"
#include "cs.h"
#include "constants.h"
#include "global_opts.h"
#include "cholmod.h"
#include <stdio.h>
#include <math.h>

csc* random_matrix(c_int m, c_int n, c_float density) {
  c_float* X_temp = c_calloc(m*n, sizeof(c_float));
  //c_int* i = c_calloc(n*m, sizeof(c_int));
  c_int* rows = c_calloc(n*m, sizeof(c_int));
  c_int* p = c_calloc(n+1, sizeof(c_int));
  int k = 0;
  int col = 0;
  c_int x_rand;
  for (int j = 0; j < n*m; j++) {
    if (j%m == 0) {
      p[col] = k;
      col++;
    }   
    x_rand = rand();
    if ((x_rand % 1000) < density * 1000) {
      X_temp[k] = (c_float) x_rand/RAND_MAX;
      rows[k] = j%m;
      k++;
    }
  }
  p[col] = k;
  c_int nnz = k;
  c_float* X = c_calloc(nnz, sizeof(c_float));
  c_int* i = c_calloc(nnz, sizeof(c_int));
  for (k = 0; k < nnz; k++) {
    X[k] = X_temp[k];
    i[k] = rows[k];
  }
  c_free(rows);
  c_free(X_temp);
  // //c_int nnz = n*m; 
  // c_int* p = c_calloc(n+1, sizeof(c_int));
  // for (int i = 0; i < n+1; i++) {
  //   p[i] = i*m;
  // }

  // c_int* i = c_calloc(n*m, sizeof(c_int));
  // for (int k = 0; k < n*m; k++) {
  //   i[k] = k%m;
  // }

  csc* M = csc_matrix(m, n, nnz, X, i, p);

  return M;
}

c_float* random_vector(c_int n) {
  c_float* X = c_calloc(n, sizeof(c_float));
  for (int i = 0; i < n; i++) {
    X[i] = (c_float) 10*rand()/RAND_MAX;
  }
  return X;
}

c_float* constant_vector(c_float c, c_int n) {
  c_float* X = c_calloc(n, sizeof(c_float));
  for (int i = 0; i < n; i++) {
    X[i] = c;
  }
  return X;
}
int main() {

// Load problem data
  
  //regular problem
  // c_float Q_x[4] =
  // { 4.00000000000000000000, 1.00000000000000000000, 1.00000000000000000000,
  //   2.00000000000000000000, };
  // c_int   Q_nnz  = 4;
  // c_int   Q_i[4] = { 0, 1, 0, 1, };
  // c_int   Q_p[3] = { 0, 2, 4, };
  // c_float q[2]   = { 1.00000000000000000000, 1.00000000000000000000, };
  // c_float A_x[4] =
  // { 1.0, 1.0, 1.0, 1.0};
  // c_int   A_nnz  = 4;
  // c_int   A_i[4] = { 0, 1, 0, 2, };
  // c_int   A_p[3] = { 0, 2, 4, };
  // c_float bmin[3]   =
  // { 1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000, };
  // c_float bmax[3] =
  // { 1.00000000000000000000, 0.69999999999999995559, 0.69999999999999995559, };
  // c_int n = 2;
  // c_int m = 3;

  // // primal infeasible problem
  // c_float Q[4] =
  // { 4.00000000000000000000, 1.00000000000000000000, 1.00000000000000000000,
  //   2.00000000000000000000, };
  // c_float q[2]   = { 1.00000000000000000000, 1.00000000000000000000, };
  // c_float A[6] =
  // { 1.00000000000000000000, 1.00000000000000000000, 1.00000000000000000000,
  //   -0.500000000000000000000, 0.23, 0.46 };
  // c_float bmin[3]   =   {1.0, 1.0, 1.0};
  // c_float bmax[3] =   {1.0, 1.0, 1.0};

  // // dual infeasible problem
  // c_float Q[4] = { 0, 0, 0, 0, };
  // c_float q[2]   = { 1.0, 0, };
  // c_float A[6] = { 10.0, 1.0, 1.0, 0.2, 0, 1 };
  // c_float bmin[3]   =   {-QPALM_INFTY, -QPALM_INFTY, 1.0};
  // c_float bmax[3] =   {1.0, 1.0, 1.0};

  c_int n = 3;
  c_int m = 5;



  // Problem settings
  QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));

  // Structures
  QPALMWorkspace *work; // Workspace
  QPALMData *data;      // QPALMData

  // Populate data
  data    = (QPALMData *)c_malloc(sizeof(QPALMData));
  data->n = n;
  data->m = m;
  //data->Q = csc_matrix(data->n, data->n, Q_nnz, Q_x, Q_i, Q_p);
  // data->Q = random_matrix(data->n, data->n, 1e-1);
  cholmod_dense *Q_dense, *A_dense;

  cholmod_common c;
  cholmod_start(&c);


  Q_dense = cholmod_ones(n, n, CHOLMOD_REAL, &c);
  data->Q = cholmod_dense_to_sparse(Q_dense, 1, &c);
  data->q = random_vector(data->n);
  // for (int i = 0; i < data->n; i++) {
  //   printf(" %.20f", data->q[i]);
  // }
  //data->q = q;
  //data->A = csc_matrix(data->m, data->n, A_nnz, A_x, A_i, A_p);
  // data->A = random_matrix(data->m, data->n, 1e-1);
  A_dense = cholmod_ones(m, n, CHOLMOD_REAL, &c);
  data->A = cholmod_dense_to_sparse(A_dense, 1, &c);
  data->bmin = constant_vector(-2, data->m);
  data->bmax = constant_vector(2, data->m);


  // Define Solver settings as default
  qpalm_set_default_settings(settings);

  // Setup workspace
  work = qpalm_setup(data, settings, &c);

  // Solve Problem
  qpalm_solve(work);

  printf("Solver status: ");
  printf(work->info->status);
  printf(" \n");
  printf("Iter: %d\n", (int) work->info->iter);
  printf("Iter Out: %d\n", (int) work->info->iter_out);
  // for (c_int i = 0; i < n; i++) {
  //   printf("Solution variable %.10f \n",work->x[i]);
  // }
  #ifdef PROFILING
  printf("Setup time: %f\n", work->info->setup_time);
  printf("Solve time: %f\n", work->info->solve_time);
  printf("Run time: %f\n", work->info->run_time);
  #endif
  // Clean workspace
  
  cholmod_free_dense(&Q_dense, &c);
  cholmod_free_dense(&A_dense, &c);
  cholmod_free_sparse(&data->Q, &c);
  cholmod_free_sparse(&data->A, &c);
  qpalm_cleanup(work);
  // csc_spfree(data->A);
  // csc_spfree(data->Q);
  c_free(data->q);
  c_free(data->bmin);
  c_free(data->bmax);
  c_free(data);
  c_free(settings);


  // cholmod_finish(&work->chol->c);
    return 0;

}