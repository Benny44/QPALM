#include "minunit.h"
#include "qpalm.h"
#include "lin_alg.h"
#include "global_opts.h"
#include "constants.h"

#define N 4
#define M 5
#define ANZMAX 4
#define QNZMAX 4

QPALMWorkspace *work; // Workspace
QPALMSettings *settings;
QPALMData *data;
cholmod_common *c;
cholmod_common common;
static c_float solution[N] = {2.0000000e+00, -6.3801365e+01, -3.3821109e+03, -6.0483288e+00};

void basic_qp_suite_setup(void) {
    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;

    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->c = 0;

    c = &common;
    CHOLMOD(start)(c);
    data->A = CHOLMOD(allocate_sparse)(data->m, data->n, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, c);
    data->Q = CHOLMOD(allocate_sparse)(data->n, data->n, QNZMAX , TRUE, TRUE, -1, CHOLMOD_REAL, c);
    CHOLMOD(finish)(c);
    
    c_float *Ax = data->A->x;
    c_int *Ai = data->A->i;
    c_int *Ap = data->A->p;
    Ax[0] = -1.0000000000000000;
    Ai[0] = 3;
    Ax[1] = 0.0254311360000000;
    Ai[1] = 4;
    Ax[2] = -0.0001000000000000;
    Ai[2] = 0;
    Ax[3] = 0.3306698500000000;
    Ai[3] = 2;
    Ap[0] = 0;
    Ap[1] = 1;
    Ap[2] = 2;
    Ap[3] = 3;
    Ap[4] = 4;

    c_float *Qx = data->Q->x;
    c_int *Qi = data->Q->i;
    c_int *Qp = data->Q->p;
    Qx[0] = 1.0000000000000000;
    Qi[0] = 0;
    Qx[1] = 0.0464158880000000;
    Qi[1] = 1;
    Qx[2] = 0.0021544347000000;
    Qi[2] = 2;
    Qx[3] = 0.0001000000000000;
    Qi[3] = 3;
    Qp[0] = 0;
    Qp[1] = 1;
    Qp[2] = 2;
    Qp[3] = 3;
    Qp[4] = 4;

    data->q = c_calloc(N,sizeof(c_float));
    data->bmin = c_calloc(M,sizeof(c_float));
    data->bmax = c_calloc(M,sizeof(c_float));
    data->bmin[0] = -2; 
    data->bmin[1] = -2;
    data->bmin[2] = -2;
    data->bmin[3] = -2;
    data->bmin[4] = -2;

    data->bmax[0] = 2; 
    data->bmax[1] = 2;
    data->bmax[2] = 2;
    data->bmax[3] = 2;
    data->bmax[4] = 2;

    data->q[0] = -2.0146781e+00; 
    data->q[1] = 2.9613971e+00; 
    data->q[2] = 7.2865370e+00; 
    data->q[3] = 7.8925204e+00;
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
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

MU_TEST(test_basic_qp_unscaled) {
    // Setup workspace
    settings->scaling = 0;
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}
MU_TEST(test_basic_qp_noprox) {
    // Setup workspace
    settings->proximal = FALSE;
    settings->scaling = 2;
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}
MU_TEST(test_basic_qp_noprox_unscaled) {
    // Setup workspace
    settings->proximal = FALSE;
    settings->scaling = 0;
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

MU_TEST(test_basic_qp_warm_start) {
    // Setup workspace
    settings->scaling = 2;
    settings->proximal = TRUE;
    settings->warm_start = TRUE;
    work = qpalm_setup(data, settings, c);
    c_float x[N] = {2.0, -60.0, -3380.0, -6.0};
    c_float y[M] = {0.0, 0.0, -23.0, -0.014, 0.0};
    qpalm_warm_start(work, x, y);

    // Solve Problem
    qpalm_solve(work);
    mu_assert_true(work->info->iter < 12);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

MU_TEST(test_basic_qp_warm_start_unscaled) {
    // Setup workspace
    settings->scaling = 0;
    settings->proximal = TRUE;
    settings->warm_start = TRUE;
    work = qpalm_setup(data, settings, c);
    c_float x[N] = {2.0, -60.0, -3380.0, -6.0};
    c_float y[M] = {0.0, 0.0, -23.0, -0.01, 0.0};
    qpalm_warm_start(work, x, y);

    // Solve Problem
    qpalm_solve(work);
    mu_assert_true(work->info->iter < 12);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

MU_TEST(test_basic_qp_warm_start_noprox) {
    // Setup workspace
    settings->scaling = 2;
    settings->proximal = FALSE;
    settings->warm_start = TRUE;
    work = qpalm_setup(data, settings, c);
    c_float x[N] = {2.0, -60.0, -3380.0, -6.0};
    c_float y[M] = {0.0, 0.0, -23.0, -0.01, 0.0};
    qpalm_warm_start(work, x, y);

    // Solve Problem
    qpalm_solve(work);
    mu_assert_true(work->info->iter < 12);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

MU_TEST(test_basic_qp_warm_start_noprox_unscaled) {
    // Setup workspace
    settings->scaling = 0;
    settings->proximal = FALSE;
    settings->warm_start = TRUE;
    work = qpalm_setup(data, settings, c);
    c_float x[N] = {2.0, -60.0, -3380.0, -6.0};
    c_float y[M] = {0.0, 0.0, -23.0, -0.01, 0.0};
    qpalm_warm_start(work, x, y);

    // Solve Problem
    qpalm_solve(work);
    mu_assert_true(work->info->iter < 12);
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
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
    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);

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
    c_float tol = 1e-15;
    for(c_int i = 0; i < N; i++) {
        mu_assert_double_eq(work->solution->x[i], x_sol[i], tol);
    }
    for(c_int i = 0; i < M; i++) {
        mu_assert_double_eq(work->solution->y[i], y_sol[i], tol);
    }

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
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
    mu_assert_double_eq(work->info->objective, work->info->dual_objective, 1e-5);
}

MU_TEST(test_basic_qp_dual_early_termination) {
    // Setup workspace
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->enable_dual_termination = TRUE;
    settings->dual_objective_limit = -1000000000.0;

    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_DUAL_TERMINATED);
    mu_assert_long_eq(work->info->iter_out, 0); /* terminate on the first outer iteration */

    settings->enable_dual_termination = FALSE;
}

MU_TEST(test_basic_qp_sigma_max) {
    // Setup workspace
    settings->sigma_max = 1e3;
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
}

MU_TEST(test_basic_qp_time_limit) {
    // Setup workspace
    settings->time_limit = 0.01*1e-3;
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_TIME_LIMIT_REACHED);
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
    MU_RUN_TEST(test_basic_qp_sigma_max);
    MU_RUN_TEST(test_basic_qp_time_limit);
}
