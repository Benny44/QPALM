#include "minunit.h"
#include "qpalm.h"
#include "lin_alg.h"
#include "global_opts.h"
#include "constants.h"

#define N 2
#define M 3
#define ANZMAX 4
#define QNZMAX 2

QPALMWorkspace *work; // Workspace
QPALMSettings *settings;
QPALMData *data;
cholmod_common *c;
cholmod_common common;

void basic_qp_suite_setup(void) {
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
}

void basic_qp_suite_teardown(void) {
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

void basic_qp_test_teardown(void) {
    qpalm_cleanup(work);
}



MU_TEST(test_basic_qp) {
    // Setup workspace
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], -0.1, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.3, 1e-5);
}

MU_TEST(test_basic_qp_unscaled) {
    // Setup workspace
    settings->scaling = 0;
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], -0.1, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.3, 1e-5);
}
MU_TEST(test_basic_qp_noprox) {
    // Setup workspace
    settings->proximal = FALSE;
    settings->scaling = 2;
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], -0.1, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.3, 1e-5);
}
MU_TEST(test_basic_qp_noprox_unscaled) {
    // Setup workspace
    settings->proximal = FALSE;
    settings->scaling = 0;
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], -0.1, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.3, 1e-5);
}

MU_TEST(test_basic_qp_warm_start) {
    // Setup workspace
    settings->scaling = 2;
    settings->proximal = TRUE;
    settings->warm_start = TRUE;
    work = qpalm_setup(data, settings, c);
    c_float x[N] = {-0.1, 0.3};
    c_float y[M] = {-0.9, 1.55, 0.0};
    qpalm_warm_start(work, x, y);

    // Solve Problem
    qpalm_solve(work);
    mu_assert_true(work->info->iter < 2);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], -0.1, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.3, 1e-5);
}

MU_TEST(test_basic_qp_warm_start_unscaled) {
    // Setup workspace
    settings->scaling = 0;
    settings->proximal = TRUE;
    settings->warm_start = TRUE;
    work = qpalm_setup(data, settings, c);
    c_float x[N] = {-0.1, 0.3};
    c_float y[M] = {-0.9, 1.55, 0.0};
    qpalm_warm_start(work, x, y);

    // Solve Problem
    qpalm_solve(work);
    mu_assert_true(work->info->iter < 2);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], -0.1, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.3, 1e-5);
}

MU_TEST(test_basic_qp_warm_start_noprox) {
    // Setup workspace
    settings->scaling = 2;
    settings->proximal = FALSE;
    settings->warm_start = TRUE;
    work = qpalm_setup(data, settings, c);
    c_float x[N] = {-0.1, 0.3};
    c_float y[M] = {-0.9, 1.55, 0.0};
    qpalm_warm_start(work, x, y);

    // Solve Problem
    qpalm_solve(work);
    mu_assert_true(work->info->iter < 2);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], -0.1, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.3, 1e-5);
}

MU_TEST(test_basic_qp_warm_start_noprox_unscaled) {
    // Setup workspace
    settings->scaling = 0;
    settings->proximal = FALSE;
    settings->warm_start = TRUE;
    work = qpalm_setup(data, settings, c);
    c_float x[N] = {-0.1, 0.3};
    c_float y[M] = {-0.9, 1.55, 0.0};
    qpalm_warm_start(work, x, y);

    // Solve Problem
    qpalm_solve(work);
    mu_assert_true(work->info->iter < 2);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], -0.1, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.3, 1e-5);
}

MU_TEST(test_basic_qp_warm_start_resolve) {
    // Setup workspace
    work = qpalm_setup(data, settings, c);
    //Store initial guesses    
    c_float x[N]; 
    c_float y[M];
    prea_vec_copy(work->x, x, N);
    prea_vec_copy(work->y, y, M);
    
    // Solve Problem
    qpalm_solve(work);
    // mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);

    // Store solution
    c_float x_sol[N], y_sol[M];
    prea_vec_copy(work->solution->x, x_sol, N);
    prea_vec_copy(work->solution->y, y_sol, M);
    c_int iter = work->info->iter;
    
    // Warm start and resolve problem
    qpalm_warm_start(work, x, y);
    qpalm_solve(work);
    // mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_long_eq(work->info->iter, iter);
    mu_assert_double_eq(work->solution->x[0], x_sol[0], 1e-15);
    mu_assert_double_eq(work->solution->x[1], x_sol[1], 1e-15);
    mu_assert_double_eq(work->solution->y[0], y_sol[0], 1e-15);
    mu_assert_double_eq(work->solution->y[1], y_sol[1], 1e-15);
    mu_assert_double_eq(work->solution->y[2], y_sol[2], 1e-15);

}

MU_TEST(test_basic_qp_maxiter) {
    // Setup workspace
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->max_iter = 1;

    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_MAX_ITER_REACHED);
}

MU_TEST(test_basic_qp_inner_maxiter) {
    // Setup workspace
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-8;
    settings->eps_rel = 1e-8;
    settings->inner_max_iter = 2;
    settings->max_iter = 10;

    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_MAX_ITER_REACHED);
}

MU_TEST(test_basic_qp_dual_objective) {
    // Setup workspace
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->enable_dual_termination = TRUE;

    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->solution->x[0], -0.1, 1e-5);
    mu_assert_double_eq(work->solution->x[1], 0.3, 1e-5);
    mu_assert_double_eq(work->info->objective, work->info->dual_objective, 1e-5);
}

MU_TEST(test_basic_qp_dual_early_termination) {
    // Setup workspace
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->enable_dual_termination = TRUE;
    settings->dual_objective_limit = -1000.0;

    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_DUAL_TERMINATED);
    mu_assert_long_eq(work->info->iter_out, 0); /* terminate on the first outer iteration */
}




MU_TEST_SUITE(suite_basic_qp) {
    MU_SUITE_CONFIGURE(basic_qp_suite_setup, basic_qp_suite_teardown, NULL, basic_qp_test_teardown);

    MU_RUN_TEST(test_basic_qp);
    MU_RUN_TEST(test_basic_qp_unscaled);
    MU_RUN_TEST(test_basic_qp_noprox);
    MU_RUN_TEST(test_basic_qp_noprox_unscaled);
    MU_RUN_TEST(test_basic_qp_warm_start);
    MU_RUN_TEST(test_basic_qp_warm_start_unscaled);
    MU_RUN_TEST(test_basic_qp_warm_start_noprox);
    MU_RUN_TEST(test_basic_qp_warm_start_noprox_unscaled);
    MU_RUN_TEST(test_basic_qp_warm_start_resolve);
    MU_RUN_TEST(test_basic_qp_maxiter);
    MU_RUN_TEST(test_basic_qp_inner_maxiter);
    MU_RUN_TEST(test_basic_qp_dual_objective);
    MU_RUN_TEST(test_basic_qp_dual_early_termination);
}
