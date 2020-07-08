#include "minunit.h"
#include "qpalm.h"

#define N 2
#define M 3
#define ANZMAX 4
#define QNZMAX 2

QPALMWorkspace *work; // Workspace
QPALMSettings *settings;
QPALMData *data;
solver_common *c;
solver_common common;


void update_suite_setup(void) {
    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->scaling = 2;
    settings->proximal = TRUE;
    settings->factorization_method = FACTORIZE_KKT;

    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->c = 0;
    data->q = c_calloc(N,sizeof(c_float));
    data->q[0] = 1; data->q[1] = -2; 
    data->bmin = c_calloc(M,sizeof(c_float));
    data->bmin[0] = -1; data->bmin[1] = -3; data->bmin[2] = -0.2; 
    data->bmax = c_calloc(M,sizeof(c_float));
    data->bmax[0] = 1; data->bmax[1] = 3; data->bmax[2] = 0.2; 

    // solver_common common;
    c = &common;
    #ifdef USE_LADEL
    solver_sparse *A = ladel_sparse_alloc(M, N, ANZMAX, UNSYMMETRIC, TRUE, FALSE);
    solver_sparse *Q = ladel_sparse_alloc(N, N, QNZMAX, UPPER, TRUE, FALSE);
    #elif defined USE_CHOLMOD
    CHOLMOD(start)(c);
    solver_sparse *A = CHOLMOD(allocate_sparse)(M, N, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, c);
    solver_sparse *Q = CHOLMOD(allocate_sparse)(N, N, QNZMAX, TRUE, TRUE, -1, CHOLMOD_REAL, c);
    CHOLMOD(finish)(c);
    #endif /* USE_CHOLMOD */
    c_float *Ax;
    c_int *Ai, *Ap;
    Ax = A->x;
    Ap = A->p;
    Ai = A->i;
    Ax[0] = 10.0; Ax[1] = 1.0; Ax[2] = 10.0; Ax[3] = 1.0;
    Ap[0] = 0; Ap[1] = 2; Ap[2] = 4;
    Ai[0] = 0; Ai[1] = 2; Ai[2] = 1; Ai[3] = 2;
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

    // Setup workspace
    work = qpalm_setup(data, settings);
}

void update_suite_teardown(void) {
    c_free(settings);
    // Clean setup
    #ifdef USE_LADEL
    data->Q = ladel_sparse_free(data->Q);
    data->A = ladel_sparse_free(data->A);
    #elif defined USE_CHOLMOD
    CHOLMOD(start)(c);
    CHOLMOD(free_sparse)(&data->Q, c);
    CHOLMOD(free_sparse)(&data->A, c);
    CHOLMOD(finish)(c);
    #endif /* USE_CHOLMOD */
    c_free(data->q);
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);

    qpalm_cleanup(work);
}

MU_TEST(test_update_settings) {
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], -0.1, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.3, 1e-5);

    settings->gamma_init *= 0.1;
    settings->theta = 0.9;
    settings->proximal = TRUE;
    settings->scaling = 10;

    qpalm_update_settings(work, settings);
    mu_assert_false(work->info->status_val == QPALM_ERROR);

    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], -0.1, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.3, 1e-5);
}

MU_TEST(test_update_bounds) {
    data->bmin[0] = 0.0;
    data->bmax[1] = 1.5;
    qpalm_update_bounds(work, data->bmin, data->bmax);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], 0.0, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.15, 1e-5);

    //reset bounds
    data->bmin[0] = -1;
    data->bmax[1] = 3;
    qpalm_update_bounds(work, data->bmin, data->bmax);
}
MU_TEST(test_update_q) {
    
    data->q[0] = -0.5; data->q[1] = -0.75;
    qpalm_update_q(work, data->q);

    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], 0.02, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.18, 1e-5);
}

MU_TEST_SUITE(suite_update){
    MU_SUITE_CONFIGURE(update_suite_setup, update_suite_teardown, NULL, NULL);

    MU_RUN_TEST(test_update_settings);
    MU_RUN_TEST(test_update_bounds);
    MU_RUN_TEST(test_update_q);
    
}
