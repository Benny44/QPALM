#include "qpalm.h"
#include "global_opts.h"
#include "constants.h"
#include <CUnit/CUnit.h>

#define N 2
#define M 3
#define ANZMAX 6
#define QNZMAX 2

// Structures
QPALMWorkspace *work; // Workspace

void dua_inf_qp_setup(void) {
    QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-8;
    settings->eps_rel = 1e-8;

    QPALMData *data    = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    c_float q[N] = {1, -2};
    data->q = q;
    c_float bmin[M] = {-5, -10, -20};
    c_float bmax[M] = {5, 10, 20};
    data->bmin = bmin;
    data->bmax = bmax;

    cholmod_common c;
    CHOLMOD(start)(&c);
    cholmod_sparse *A = CHOLMOD(allocate_sparse)(M, N, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    c_float *Ax;
    c_int *Ai, *Ap;
    Ax = A->x;
    Ap = A->p;
    Ai = A->i;
    Ax[0] = 1.0; Ax[1] = 1.0; Ax[2] = 1.0; Ax[3] = 1.0; Ax[4] = 1.0; Ax[5] = 1.0; 
    Ap[0] = 0; Ap[1] = 3; Ap[2] = 6;
    Ai[0] = 0; Ai[1] = 1; Ai[2] = 2; Ai[3] = 0; Ai[4] = 1; Ai[5] = 2;


    cholmod_sparse *Q = CHOLMOD(allocate_sparse)(N, N, QNZMAX, TRUE, TRUE, -1, CHOLMOD_REAL, &c);
    c_float *Qx;
    c_int *Qi, *Qp;
    Qx = Q->x;
    Qp = Q->p;
    Qi = Q->i;
    Qx[0] = 0.0; Qx[1] = 0.0; 
    Qp[0] = 0; Qp[1] = 1; Qp[2] = 2;
    Qi[0] = 0; Qi[1] = 1; 

    data->A = A;
    data->Q = Q;
    CHOLMOD(finish)(&c);
    // Setup workspace
    work = qpalm_setup(data, settings, &c);

    // Clean setup
    CHOLMOD(start)(&c);
    CHOLMOD(free_sparse)(&data->Q, &c);
    CHOLMOD(free_sparse)(&data->A, &c);
    CHOLMOD(finish)(&c);

    c_free(data);
    c_free(settings);
}

void dua_inf_qp_teardown(void) {
    qpalm_cleanup(work);
}

void test_dua_inf_qp(void) {
    // Solve Problem
    qpalm_solve(work);

    CU_ASSERT_EQUAL(work->info->status_val, QPALM_DUAL_INFEASIBLE);
}
