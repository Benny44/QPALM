#include "minunit.h"
#include "qpalm.h"
#include "lin_alg.h"
#include "global_opts.h"
#include "constants.h"
#include "util.h"
#ifdef USE_LADEL
#include "ladel_global.h"
#elif defined USE_CHOLMOD
#include "cholmod.h"
#endif

#define N 2
#define M 2
#define ANZMAX 2
#define QNZMAX 2

QPALMWorkspace *work; // Workspace
QPALMSettings *settings;
QPALMData *data;
solver_common *c;
solver_common common;

/* Problem to exercise the linesearch when all breakpoints are traversed. */

void ls_qp_suite_setup(void) {
    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->gamma_max = 1e3;

    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->c = 0;

    c = &common;
    #ifdef USE_LADEL
    data->A = ladel_sparse_alloc(M, N, ANZMAX, UNSYMMETRIC, TRUE);
    data->Q = ladel_sparse_alloc(N, N, QNZMAX, UPPER, TRUE);
    #elif defined USE_CHOLMOD
    CHOLMOD(start)(c);
    data->A = CHOLMOD(allocate_sparse)(M, N, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, c);
    data->Q = CHOLMOD(allocate_sparse)(N, N, QNZMAX, TRUE, TRUE, -1, CHOLMOD_REAL, c);
    CHOLMOD(finish)(c);
    #endif /* USE_CHOLMOD */

    c_float *Ax = data->A->x;
    c_int *Ai = data->A->i;
    c_int *Ap = data->A->p;
    Ax[0] = -1.0000000000000000;
    Ai[0] = 1;
    Ax[1] = 0.0001000000000000;
    Ai[1] = 0;
    Ap[0] = 0;
    Ap[1] = 1;
    Ap[2] = 2;

    c_float *Qx = data->Q->x;
    c_int *Qi = data->Q->i;
    c_int *Qp = data->Q->p;
    Qx[0] = 1.0000000000000000;
    Qi[0] = 0;
    Qx[1] = 0.0001000000000000;
    Qi[1] = 1;
    Qp[0] = 0;
    Qp[1] = 1;
    Qp[2] = 2;

    data->q = c_calloc(N,sizeof(c_float));
    data->bmin = c_calloc(M,sizeof(c_float));
    data->bmax = c_calloc(M,sizeof(c_float));
    data->bmin[0] = -2; data->bmin[1] = -2;
    data->bmax[0] = 2; data->bmax[1] = 2;
    data->q[0] = 2.5150105e+00; data->q[1] = 1.6259589e+01;
}

void ls_qp_suite_teardown(void) {
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

void ls_qp_test_teardown(void) {
    qpalm_cleanup(work);
}


MU_TEST(test_ls_qp) {
    // Setup workspace
    work = qpalm_setup(data, settings);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);

    c_float *solution = c_calloc(N,sizeof(c_float));
    solution[0] = -2.0000000e+00;
    solution[1] = -2.0000000e+04;
    

    for(c_int i = 0; i < N; i++) {
        mu_assert_double_eq(work->solution->x[i], solution[i], 1e-5);
    }

    c_free(solution);
}

MU_TEST_SUITE(suite_ls_qp) {
    MU_SUITE_CONFIGURE(ls_qp_suite_setup, ls_qp_suite_teardown, NULL, ls_qp_test_teardown);

    MU_RUN_TEST(test_ls_qp);
}
