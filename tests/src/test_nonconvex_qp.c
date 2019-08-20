#include "minunit.h"
#include "qpalm.h"
#include "global_opts.h"
#include "constants.h"

#define N 2
#define M 2
#define ANZMAX 2
#define QNZMAX 2

QPALMWorkspace *work; // Workspace

void nonconvex_qp_suite_setup(void) {
    QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->nonconvex = TRUE;

    QPALMData *data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->q = c_calloc(N,sizeof(c_float));
    data->q[0] = 1; data->q[1] = 2; 
    data->bmin = c_calloc(M,sizeof(c_float));
    data->bmin[0] = -2; data->bmin[1] = -2; 
    data->bmax = c_calloc(M,sizeof(c_float));
    data->bmax[0] = 2; data->bmax[1] = 2;

    cholmod_common common, *c;
    c = &common;
    CHOLMOD(start)(c);
    cholmod_sparse *A = CHOLMOD(allocate_sparse)(M, N, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, c);
    c_float *Ax;
    c_int *Ai, *Ap;
    Ax = A->x;
    Ap = A->p;
    Ai = A->i;
    Ax[0] = 1.0; Ax[1] = 1.0; 
    Ap[0] = 0; Ap[1] = 1; Ap[2] = 2;
    Ai[0] = 0; Ai[1] = 1; 
    cholmod_sparse *Q = CHOLMOD(allocate_sparse)(N, N, QNZMAX, TRUE, TRUE, -1, CHOLMOD_REAL, c);
    c_float *Qx;
    c_int *Qi, *Qp;
    Qx = Q->x;
    Qp = Q->p;
    Qi = Q->i;
    Qx[0] = 2; Qx[1] = -3; 
    Qp[0] = 0; Qp[1] = 1; Qp[2] = 2;
    Qi[0] = 0; Qi[1] = 1; 

    data->A = A;
    data->Q = Q;
    CHOLMOD(finish)(c); 

    // Setup workspace
    work = qpalm_setup(data, settings, c);

    // Cleanup temporary structures
    c_free(settings);
    CHOLMOD(start)(c);
    CHOLMOD(free_sparse)(&data->Q, c);
    CHOLMOD(free_sparse)(&data->A, c);
    CHOLMOD(finish)(c);
    c_free(data->q);
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);
}

void nonconvex_qp_suite_teardown(void) {
    qpalm_cleanup(work);
}

void test_nonconvex_qp(void) {
    // Solve Problem
    qpalm_solve(work);

    mu_assert_cint_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], -0.5, 1e-5);
    mu_assert_double_eq(work->solution->x[1], -2.0, 1e-5);
}

MU_TEST_SUITE(suite_nonconvex) {
    MU_SUITE_CONFIGURE(nonconvex_qp_suite_setup, nonconvex_qp_suite_teardown, NULL, NULL);

    MU_RUN_TEST(test_nonconvex_qp);
}
