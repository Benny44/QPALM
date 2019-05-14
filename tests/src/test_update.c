#include "qpalm.h"
#include "global_opts.h"
#include "constants.h"
#include <CUnit/CUnit.h>

#define N 2
#define M 3
#define ANZMAX 4
#define QNZMAX 2

QPALMWorkspace *work; // Workspace
QPALMSettings *settings;
QPALMData *data;
cholmod_common *c;
cholmod_common common;


int update_suite_setup(void) {
    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;

    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->q = c_calloc(N,sizeof(c_float));
    data->q[0] = 1; data->q[1] = -2; 
    data->bmin = c_calloc(M,sizeof(c_float));
    data->bmin[0] = -0.1; data->bmin[1] = -0.3; data->bmin[2] = -0.2; 
    data->bmax = c_calloc(M,sizeof(c_float));
    data->bmax[0] = 0.1; data->bmax[1] = 0.3; data->bmax[2] = 0.2; 

    // cholmod_common common;
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

    // Setup workspace
    work = qpalm_setup(data, settings, c);
    return 0;
}

int update_suite_teardown(void) {
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

    qpalm_cleanup(work);
    return 0;
}

void test_update_settings(void) {
    // Solve Problem
    qpalm_solve(work);

    CU_ASSERT_EQUAL(work->info->status_val, QPALM_SOLVED);
    CU_ASSERT_DOUBLE_EQUAL(work->solution->x[0], -0.1, 1e-5);
    CU_ASSERT_DOUBLE_EQUAL(work->solution->x[1], 0.3, 1e-5);

    settings->gamma_init *= 10;
    settings->theta = 0.9;
    settings->proximal = FALSE;

    qpalm_update_settings(work, settings);
    qpalm_solve(work);

    CU_ASSERT_EQUAL(work->info->status_val, QPALM_SOLVED);
    CU_ASSERT_DOUBLE_EQUAL(work->solution->x[0], -0.1, 1e-5);
    CU_ASSERT_DOUBLE_EQUAL(work->solution->x[1], 0.3, 1e-5);
}

void test_update_bounds(void) {
    
    data->bmin[0] = 0.0;
    data->bmax[1] = 0.15;
    qpalm_update_bounds(work, data->bmin, data->bmax);
    // Solve Problem
    qpalm_solve(work);

    CU_ASSERT_EQUAL(work->info->status_val, QPALM_SOLVED);
    CU_ASSERT_DOUBLE_EQUAL(work->solution->x[0], 0.0, 1e-5);
    CU_ASSERT_DOUBLE_EQUAL(work->solution->x[1], 0.15, 1e-5);
}
void test_update_q(void) {
    
    data->q[0] = -10;
    qpalm_update_q(work, data->q);

    // Solve Problem
    qpalm_solve(work);

    CU_ASSERT_EQUAL(work->info->status_val, QPALM_SOLVED);
    CU_ASSERT_DOUBLE_EQUAL(work->solution->x[0], 0.1, 1e-5);
    CU_ASSERT_DOUBLE_EQUAL(work->solution->x[1], 0.1, 1e-5);
}
