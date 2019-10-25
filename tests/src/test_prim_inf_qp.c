#include "minunit.h"
#include "qpalm.h"
#include "global_opts.h"
#include "constants.h"

#define N 2
#define M 3
#define ANZMAX 4
#define QNZMAX 2

// Structures
QPALMWorkspace *work; // Workspace
QPALMSettings *settings;
QPALMData *data;
cholmod_common *c;
cholmod_common common;


void prim_inf_qp_suite_setup(void) {
    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;

    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->c = 0;
    // c_float q[N] = {1, -2};
    // data->q = q;
    // c_float bmin[M] = {-5, -10, 16};
    // c_float bmax[M] = {5, 10, 20};
    // data->bmin = bmin;
    // data->bmax = bmax;

    data->q = c_calloc(N,sizeof(c_float));
    data->q[0] = 1; data->q[1] = -2; 
    data->bmin = c_calloc(M,sizeof(c_float));
    data->bmin[0] = -5; data->bmin[1] = -10; data->bmin[2] = 16; 
    data->bmax = c_calloc(M,sizeof(c_float));
    data->bmax[0] = 5; data->bmax[1] = 10; data->bmax[2] = 20; 

    c = &common;
    CHOLMOD(start)(c);
    cholmod_sparse *A = CHOLMOD(allocate_sparse)(M, N, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, c);
    c_float *Ax;
    c_int *Ai, *Ap;
    Ax = A->x;
    Ap = A->p;
    Ai = A->i;
    Ax[0] = 1.0; Ax[1] = 1.0; Ax[2] = 1.0; Ax[3] = 1.0; 
    Ap[0] = 0; Ap[1] = 2; Ap[2] = 4;
    Ai[0] = 0; Ai[1] = 2; Ai[2] = 1; Ai[3] = 2;

    cholmod_sparse *Q = CHOLMOD(allocate_sparse)(N, N, QNZMAX, TRUE, TRUE, -1, CHOLMOD_REAL, c);
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
    CHOLMOD(finish)(c);

}

void prim_inf_qp_suite_teardown(void) {
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


void prim_inf_qp_test_teardown(void) {
    qpalm_cleanup(work);
}

MU_TEST(test_prim_inf_qp) {
    settings->proximal = TRUE;
    settings->scaling = 2;
    // Setup workspace
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_PRIMAL_INFEASIBLE);
}
MU_TEST(test_prim_inf_qp_unscaled) {
    settings->proximal = TRUE;
    settings->scaling = 0;
    // Setup workspace
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_PRIMAL_INFEASIBLE);
}
MU_TEST(test_prim_inf_qp_noprox) {
    settings->proximal = FALSE;
    settings->scaling = 2;
    // Setup workspace
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_PRIMAL_INFEASIBLE);
}
MU_TEST(test_prim_inf_qp_noprox_unscaled) {
    settings->proximal = FALSE;
    settings->scaling = 0;
    // Setup workspace
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_PRIMAL_INFEASIBLE);
}

MU_TEST_SUITE(suite_prim_inf_qp) {
    MU_SUITE_CONFIGURE(prim_inf_qp_suite_setup, prim_inf_qp_suite_teardown, NULL, prim_inf_qp_test_teardown);

    MU_RUN_TEST(test_prim_inf_qp);
    MU_RUN_TEST(test_prim_inf_qp_unscaled);
    MU_RUN_TEST(test_prim_inf_qp_noprox);
    MU_RUN_TEST(test_prim_inf_qp_noprox_unscaled);
    
}
