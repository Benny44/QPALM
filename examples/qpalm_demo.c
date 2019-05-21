#include "qpalm.h"
#include "constants.h"
#include "global_opts.h"
#include "cholmod.h"
#include <stdio.h>
#include <math.h>

#define N 2
#define M 3
#define ANZMAX 4
#define QNZMAX 2

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

  c_int n = N;
  c_int m = M;

  // Problem settings
  QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));

  // Structures
  QPALMWorkspace *work; // Workspace
  QPALMData *data;      // QPALMData

  // Populate data
  data    = (QPALMData *)c_malloc(sizeof(QPALMData));
  data->n = n;
  data->m = m;
  data->q = random_vector(data->n);
  data->bmin = constant_vector(-2, data->m);
  data->bmax = constant_vector(2, data->m);

  cholmod_common c;
  CHOLMOD(start)(&c);
  cholmod_sparse *A = CHOLMOD(allocate_sparse)(m, n, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
  c_float *Ax;
  c_int *Ai, *Ap;
  Ax = A->x;
  Ap = A->p;
  Ai = A->i;
  Ax[0] = 1.0; Ax[1] = 1.0; Ax[2] = 1.0; Ax[3] = 1.0;
  Ap[0] = 0; Ap[1] = 2; Ap[2] = 4;
  Ai[0] = 0; Ai[1] = 2; Ai[2] = 1; Ai[3] = 2;
  cholmod_sparse *Q = CHOLMOD(allocate_sparse)(N, N, QNZMAX, TRUE, TRUE, -1, CHOLMOD_REAL, &c);
  c_float *Qx;
  c_int *Qi, *Qp;
  Qx = Q->x;
  Qp = Q->p;
  Qi = Q->i;
  Qx[0] = 1.0; Qx[1] = 1.5; 
  Qp[0] = 0; Qp[1] = 1; Qp[2] = 2;
  Qi[0] = 0; Qi[1] = 1; 

  data->A = A;
  data->Q = Q;

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