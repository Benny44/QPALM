#include "minunit.h"
#include "qpalm.h"
#include "lin_alg.h"
#include "global_opts.h"
#include "constants.h"
#include "util.h"
#ifdef USE_LADEL
#include "ladel_global.h"
#endif

#define N 2
#define M 3
#define ANZMAX 4
#define QNZMAX 2

QPALMWorkspace *work; // Workspace
QPALMSettings *settings;
QPALMData *data;
solver_common *c;
solver_common common;

void error_handling_suite_setup(void) {
    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);

    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->c = 0;
    data->q = c_calloc(N,sizeof(c_float));
    data->q[0] = 1; data->q[1] = -2; 
    data->bmin = c_calloc(M,sizeof(c_float));
    data->bmin[0] = -0.1; data->bmin[1] = -0.3; data->bmin[2] = -0.2; 
    data->bmax = c_calloc(M,sizeof(c_float));
    data->bmax[0] = 0.1; data->bmax[1] = 0.3; data->bmax[2] = 0.2; 

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
    Ax[0] = 1.0; Ax[1] = 1.0; Ax[2] = 1.0; Ax[3] = 1.0;
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
}

void error_handling_suite_teardown(void) {
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

void error_handling_test_teardown(void) {
    qpalm_cleanup(work);
}


MU_TEST(test_invalid_settings_during_setup) {
    // Setup workspace
    settings->max_iter = -1;
    work = qpalm_setup(data, settings);
    mu_assert_true(work == QPALM_NULL);

    settings->max_iter = 1000;
}

MU_TEST(test_invalid_data_during_setup) {
    // Setup workspace
    data->bmin[0] = 5; data->bmax[0] = 0;
    work = qpalm_setup(data, settings);
    mu_assert_true(work == QPALM_NULL);

    data->bmin[0] = -0.1;
}

MU_TEST(test_invalid_settings_during_update_settings) {
    work = qpalm_setup(data, settings);
    mu_assert_long_eq(work->info->status_val, QPALM_UNSOLVED);

    settings->max_iter = -10;
    qpalm_update_settings(work, settings);
    mu_assert_long_eq(work->info->status_val, QPALM_ERROR);

    settings->max_iter = 1000;
}

MU_TEST(test_invalid_scaling_during_update_settings) {
    work = qpalm_setup(data, settings);
    mu_assert_long_eq(work->info->status_val, QPALM_UNSOLVED);

    /*ask to unscale during update -> not possible to decrease scaling iterations */
    settings->scaling = 0;
    qpalm_update_settings(work, settings);
    mu_assert_long_eq(work->info->status_val, QPALM_ERROR);

    settings->scaling = 2;
}

MU_TEST(test_invalid_data_during_update_bounds) {
    work = qpalm_setup(data, settings);
    mu_assert_long_eq(work->info->status_val, QPALM_UNSOLVED);

    data->bmin[0] = 5; data->bmax[0] = 0;
    qpalm_update_bounds(work, data->bmin, data->bmax);
    mu_assert_long_eq(work->info->status_val, QPALM_ERROR);

    data->bmin[0] = -5;
}

MU_TEST(test_invalid_status_value) {
    work = qpalm_setup(data, settings);
    update_status(work->info, 999);
    mu_assert_string_eq(work->info->status, "unrecognised status value");

    update_status(work->info, QPALM_UNSOLVED); /*reset the status*/

    print_final_message(work);
    mu_assert_string_eq(work->info->status, "unrecognised status value");
}


MU_TEST_SUITE(suite_error_handling) {
    MU_SUITE_CONFIGURE(error_handling_suite_setup, error_handling_suite_teardown, NULL, error_handling_test_teardown);

    MU_RUN_TEST(test_invalid_settings_during_setup);
    MU_RUN_TEST(test_invalid_data_during_setup);
    MU_RUN_TEST(test_invalid_settings_during_update_settings);
    MU_RUN_TEST(test_invalid_scaling_during_update_settings);
    MU_RUN_TEST(test_invalid_data_during_update_bounds);
    MU_RUN_TEST(test_invalid_status_value);
}
