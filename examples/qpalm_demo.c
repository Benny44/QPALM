#include "qpalm.h"
#include "constants.h"
#include "global_opts.h"
#include "cholmod.h"
#include <stdio.h>
#include <math.h>

// csc* random_matrix(c_int m, c_int n, c_float density) {
//   c_float* X_temp = c_calloc(m*n, sizeof(c_float));
//   //c_int* i = c_calloc(n*m, sizeof(c_int));
//   c_int* rows = c_calloc(n*m, sizeof(c_int));
//   c_int* p = c_calloc(n+1, sizeof(c_int));
//   int k = 0;
//   int col = 0;
//   c_int x_rand;
//   for (int j = 0; j < n*m; j++) {
//     if (j%m == 0) {
//       p[col] = k;
//       col++;
//     }   
//     x_rand = rand();
//     if ((x_rand % 1000) < density * 1000) {
//       X_temp[k] = (c_float) x_rand/RAND_MAX;
//       rows[k] = j%m;
//       k++;
//     }
//   }
//   p[col] = k;
//   c_int nnz = k;
//   c_float* X = c_calloc(nnz, sizeof(c_float));
//   c_int* i = c_calloc(nnz, sizeof(c_int));
//   for (k = 0; k < nnz; k++) {
//     X[k] = X_temp[k];
//     i[k] = rows[k];
//   }
//   c_free(rows);
//   c_free(X_temp);
//   csc* M = csc_matrix(m, n, nnz, X, i, p);

//   return M;
// }

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

  c_int n = 20;
  c_int m = 50;

  // Problem settings
  QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));

  // Structures
  QPALMWorkspace *work; // Workspace
  QPALMData *data;      // QPALMData

  // Populate data
  data    = (QPALMData *)c_malloc(sizeof(QPALMData));
  data->n = n;
  data->m = m;

  cholmod_dense *Q_dense, *A_dense;

  cholmod_common c;
  CHOLMOD(start)(&c);


  Q_dense = CHOLMOD(ones)(n, n, CHOLMOD_REAL, &c);
  data->Q = CHOLMOD(dense_to_sparse)(Q_dense, 1, &c);
  data->q = random_vector(data->n);

  A_dense = CHOLMOD(ones)(m, n, CHOLMOD_REAL, &c);
  data->A = CHOLMOD(dense_to_sparse)(A_dense, 1, &c);
  data->bmin = constant_vector(-2, data->m);
  data->bmax = constant_vector(2, data->m);

  CHOLMOD(finish)(&c);

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
  CHOLMOD(start)(&work->chol->c);
  CHOLMOD(free_dense)(&Q_dense, &c);
  CHOLMOD(free_dense)(&A_dense, &c);
  CHOLMOD(free_sparse)(&data->Q, &c);
  CHOLMOD(free_sparse)(&data->A, &c);
  CHOLMOD(finish)(&work->chol->c);

  qpalm_cleanup(work);
  c_free(data->q);
  c_free(data->bmin);
  c_free(data->bmax);
  c_free(data);
  c_free(settings);

    return 0;

}