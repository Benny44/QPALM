#include "minunit.h"
#include "qpalm.h"
#include "lin_alg.h"
#include "global_opts.h"
#include "constants.h"
#include "util.h"
#include "cholmod.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "load_data.h"

#define N 4
#define M 5
#define ANZMAX 4
#define QNZMAX 4

QPALMWorkspace *work; // Workspace
QPALMSettings *settings;
QPALMData *data;
cholmod_common *c;
cholmod_common common;

void big_qp_suite_setup(void) {
    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->eps_dual_inf = 1e-3;

    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;

    c = &common;
    CHOLMOD(start)(c);
    data->A = CHOLMOD(allocate_sparse)(data->m, data->n, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, c);
    data->Q = CHOLMOD(allocate_sparse)(data->n, data->n, QNZMAX , TRUE, TRUE, -1, CHOLMOD_REAL, c);
    CHOLMOD(finish)(c);
    
    load_data("../../tests/data/big_qp", data);

}

void big_qp_suite_teardown(void) {
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

void big_qp_test_teardown(void) {
    qpalm_cleanup(work);
}


MU_TEST(test_big_qp) {
    // Setup workspace
    work = qpalm_setup(data, settings, c);
    // Solve Problem

    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);

    FILE* fp = fopen("../../tests/data/big_qp_solution", "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", "big_qp_solution");
    }
    c_float *solution = load_dense(fp, N);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }
    fclose(fp);
    c_free(solution);
    
}

MU_TEST(test_big_qp_noprox) {
    // Setup workspace
    settings->proximal = FALSE;
    work = qpalm_setup(data, settings, c);
    // Solve Problem

    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);

    FILE* fp = fopen("../../tests/data/big_qp_solution", "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", "big_qp_solution");
    }
    c_float *solution = load_dense(fp, N);
    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }

    fclose(fp);
    c_free(solution);
}

MU_TEST(test_big_qp_noscaling) {
    // Setup workspace
    settings->scaling = 0;
    settings->proximal = FALSE;
    work = qpalm_setup(data, settings, c);
    // Solve Problem

    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);

    FILE* fp = fopen("../../tests/data/big_qp_solution", "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", "big_qp_solution");
    }
    c_float *solution = load_dense(fp, N);

    c_float tol;
    for(c_int i = 0; i < N; i++) {
        tol = c_absval(1e-5*solution[i]);
        mu_assert_double_eq(work->solution->x[i], solution[i], tol);
    }

    fclose(fp);
    c_free(solution);
    
}

MU_TEST_SUITE(suite_big_qp) {
    MU_SUITE_CONFIGURE(big_qp_suite_setup, big_qp_suite_teardown, NULL, big_qp_test_teardown);

    MU_RUN_TEST(test_big_qp);
    MU_RUN_TEST(test_big_qp_noprox);
    MU_RUN_TEST(test_big_qp_noscaling);
}
