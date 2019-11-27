#include "minunit.h"
#include "qpalm.h"
#include "global_opts.h"
#include "constants.h"

#define N 4
#define M 5
#define ANZMAX 4
#define QNZMAX 4

QPALMWorkspace *work; // Workspace

void nonconvex_qp_suite_setup(void) {
    QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->nonconvex = TRUE;
    settings->scaling = FALSE; /*So we can retrieve the actual eigenvalue*/
 

    QPALMData *data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->c = 0;

    cholmod_common common, *c;
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
    Qx[2] = -0.0021544347000000;
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

    // Setup workspace
    work = qpalm_setup(data, settings);

    // Cleanup temporary structures
    c_free(settings);
    CHOLMOD(start)(c);
    CHOLMOD(free_sparse)(&data->Q, c);
    CHOLMOD(free_sparse)(&data->A, c);
    CHOLMOD(finish)(c);
    c_free(data->q);
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);
}

void nonconvex_qp_suite_teardown(void) {
    qpalm_cleanup(work);
}

void test_nonconvex_qp(void) {
    // Solve Problem
    qpalm_solve(work);

    mu_assert_long_eq(work->info->status_val, QPALM_SOLVED);
    mu_assert_double_eq(work->gamma, 1.0/0.0021544347000000, 1e-1*1.0/0.0021544347000000); /* inverse of the lowest eigenvalue */
    mu_assert_true(1/work->gamma > 0.0021544347000000); /* eigenvalue is underapproximated (here negative signs were dropped) */
}

MU_TEST_SUITE(suite_nonconvex) {
    MU_SUITE_CONFIGURE(nonconvex_qp_suite_setup, nonconvex_qp_suite_teardown, NULL, NULL);

    MU_RUN_TEST(test_nonconvex_qp);
}
