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

#define N 2
#define M 2
#define ANZMAX 2
#define QNZMAX 2

QPALMWorkspace *work; // Workspace
QPALMSettings *settings;
QPALMData *data;
cholmod_common *c;
cholmod_common common;

void ls_qp_suite_setup(void) {
    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->gamma_max = 1e3;

    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;

    c = &common;
    CHOLMOD(start)(c);
    data->A = CHOLMOD(allocate_sparse)(data->m, data->n, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, c);
    data->Q = CHOLMOD(allocate_sparse)(data->n, data->n, QNZMAX , TRUE, TRUE, -1, CHOLMOD_REAL, c);
    CHOLMOD(finish)(c);
    
    load_data("../../tests/data/ls_qp", data);

}

void ls_qp_suite_teardown(void) {
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

void ls_qp_test_teardown(void) {
    qpalm_cleanup(work);
}


MU_TEST(test_ls_qp) {
    // Setup workspace
    work = qpalm_setup(data, settings, c);
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);

    FILE* fp = fopen("../../tests/data/ls_qp_solution", "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", "ls_qp_solution");
    }
    c_float *solution = load_dense(fp, N);

    for(c_int i = 0; i < N; i++) {
        mu_assert_double_eq(work->solution->x[i], solution[i], 1e-5);
    }
    mu_assert_long_eq(work->info->iter, 16);

    c_free(solution);
    fclose(fp);
}

MU_TEST_SUITE(suite_ls_qp) {
    MU_SUITE_CONFIGURE(ls_qp_suite_setup, ls_qp_suite_teardown, NULL, ls_qp_test_teardown);

    MU_RUN_TEST(test_ls_qp);
}
