#include "minunit.h"
#include "qpalm.h"

#define N 2
#define M 3
#define ANZMAX 6
#define QNZMAX 2

// Structures
QPALMWorkspace *work; // Workspace
QPALMSettings *settings;
QPALMData *data;
solver_common *c;
solver_common common;

void dua_inf_qp_suite_setup(void) {
    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;

    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->c = 0;

    data->q = c_calloc(N,sizeof(c_float));
    data->q[0] = 1; data->q[1] = -2; 
    data->bmin = c_calloc(M,sizeof(c_float));
    data->bmin[0] = -5; data->bmin[1] = -10; data->bmin[2] = -20; 
    data->bmax = c_calloc(M,sizeof(c_float));
    data->bmax[0] = 5; data->bmax[1] = 10; data->bmax[2] = 20; 
    
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
    Ax[0] = 1.0; Ax[1] = 1.0; Ax[2] = 1.0; Ax[3] = 1.0; Ax[4] = 1.0; Ax[5] = 1.0; 
    Ap[0] = 0; Ap[1] = 3; Ap[2] = 6;
    Ai[0] = 0; Ai[1] = 1; Ai[2] = 2; Ai[3] = 0; Ai[4] = 1; Ai[5] = 2;


    c_float *Qx;
    c_int *Qi, *Qp;
    Qx = Q->x;
    Qp = Q->p;
    Qi = Q->i;
    Qx[0] = 1e-10; Qx[1] = 1e-10; 
    Qp[0] = 0; Qp[1] = 1; Qp[2] = 2;
    Qi[0] = 0; Qi[1] = 1; 

    data->A = A;
    data->Q = Q;
}

void dua_inf_qp_suite_teardown(void) {
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
}


void dua_inf_qp_test_teardown(void) {
    qpalm_cleanup(work);
}

MU_TEST(test_dua_inf_qp) {
    settings->proximal = TRUE;
    settings->scaling = 2;
    // Setup workspace
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_DUAL_INFEASIBLE);
}
MU_TEST(test_dua_inf_qp_unscaled) {
    settings->proximal = TRUE;
    settings->scaling = 0;
    // Setup workspace
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_DUAL_INFEASIBLE);
}
MU_TEST(test_dua_inf_qp_noprox) {
    //This will crash actually, hence the large gamma value
    // settings->proximal = FALSE;
    settings->proximal = FALSE;
    // settings->gamma_init = 1e13;
    // settings->gamma_max = settings->gamma_init;
    settings->scaling = 2;
    // Setup workspace
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_DUAL_INFEASIBLE);
}
MU_TEST(test_dua_inf_qp_noprox_unscaled) {
    //This will crash actually, hence the large gamma value
    // settings->proximal = FALSE;
    settings->proximal = FALSE;
    // settings->gamma_init = 1e13;
    // settings->gamma_max = settings->gamma_init;
    settings->scaling = 0;
    // Setup workspace
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_DUAL_INFEASIBLE);
}

MU_TEST_SUITE(suite_dua_inf_qp) {
    MU_SUITE_CONFIGURE(dua_inf_qp_suite_setup, dua_inf_qp_suite_teardown, NULL, dua_inf_qp_test_teardown);

    MU_RUN_TEST(test_dua_inf_qp);
    MU_RUN_TEST(test_dua_inf_qp_unscaled);
    MU_RUN_TEST(test_dua_inf_qp_noprox);
    MU_RUN_TEST(test_dua_inf_qp_noprox_unscaled);
    
}
