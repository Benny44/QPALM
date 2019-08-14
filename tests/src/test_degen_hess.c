#include "minunit.h"
#include "qpalm.h"
#include "global_opts.h"
#include "constants.h"

#define N 3
#define M 4
#define ANZMAX 6
#define QNZMAX 4
#define TOL 1e-5

QPALMWorkspace *work; // Workspace
QPALMSettings *settings;
QPALMData *data;
cholmod_common *c;
cholmod_common common;


void degen_hess_suite_setup(void) {
    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;

    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->q = c_calloc(N,sizeof(c_float));
    data->q[0] = -2; data->q[1] = -6; data->q[2] = 1;  
    data->bmin = c_calloc(M,sizeof(c_float));
    data->bmin[0] = 0.5; data->bmin[1] = -10; data->bmin[2] = -10; data->bmin[3] = -10; 
    data->bmax = c_calloc(M,sizeof(c_float));
    data->bmax[0] = 0.5; data->bmax[1] = 10; data->bmax[2] = 10; data->bmax[3] = 10;

    // cholmod_common common;
    c = &common;
    CHOLMOD(start)(c);
    cholmod_sparse *A = CHOLMOD(allocate_sparse)(M, N, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, c);
    c_float *Ax;
    c_int *Ai, *Ap;
    Ax = A->x;
    Ap = A->p;
    Ai = A->i;
    Ax[0] = 1.0; Ax[1] = 1.0; Ax[2] = 1.0; Ax[3] = 1.0; Ax[4] = 1.0; Ax[5] = 1.0;
    Ap[0] = 0; Ap[1] = 2; Ap[2] = 4; Ap[3] = 6;
    Ai[0] = 0; Ai[1] = 1; Ai[2] = 0; Ai[3] = 2; Ai[4] = 0; Ai[5] = 3;
    cholmod_sparse *Q = CHOLMOD(allocate_sparse)(N, N, QNZMAX, TRUE, TRUE, -1, CHOLMOD_REAL, c);
    c_float *Qx;
    c_int *Qi, *Qp;
    Qx = Q->x;
    Qp = Q->p;
    Qi = Q->i;
    Qx[0] = 1.0; Qx[1] = -1.0; Qx[2] = -1.0; Qx[3] = 2.0;  
    Qp[0] = 0; Qp[1] = 2; Qp[2] = 4; Qp[3] = 4;
    Qi[0] = 0; Qi[1] = 1; Qi[2] = 0; Qi[3] = 1; 

    data->A = A;
    data->Q = Q;
    CHOLMOD(finish)(c); 
}

void degen_hess_suite_teardown(void) {
    c_free(settings);
    // Clean setup
    CHOLMOD(start)(c);
    CHOLMOD(free_sparse)(&data->Q, c);
    CHOLMOD(free_sparse)(&data->A, c);
    CHOLMOD(finish)(c);
    c_free(data->q);
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);
}

void degen_hess_test_teardown(void) {
    qpalm_cleanup(work);
}


MU_TEST(test_degen_hess) {
    // Setup workspace
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_int_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], 5.5, TOL);
    mu_assert_double_eq(work->solution->x[1], 5, TOL);
    mu_assert_double_eq(work->solution->x[2], -10, TOL);
}

MU_TEST_SUITE(suite_degen_hess) {
    MU_SUITE_CONFIGURE(degen_hess_suite_setup, degen_hess_suite_teardown, NULL, degen_hess_test_teardown);

    MU_RUN_TEST(test_degen_hess);
}       
