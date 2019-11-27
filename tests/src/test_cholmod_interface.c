#include "minunit.h"
#include "cholmod_interface.h"
#include "global_opts.h"
#include "qpalm.h"
#include "constants.h"

#define N 2
#define M 3
#define ANZMAX 5
#define QNZMAX 4
#define TOL 1e-8

QPALMWorkspace *work; // Workspace
cholmod_sparse *A; // MxN matrix
cholmod_sparse *Q; //NxN symmetric matrix
cholmod_common common, *c;


void cholmod_suite_setup(void) {
    QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;

    QPALMData *data    = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    data->c = 0;
    c_float q[N] = {1, -2};
    data->q = q;
    c_float bmin[M] = {-5, -10, -2};
    c_float bmax[M] = {5, 10, 3};
    data->bmin = bmin;
    data->bmax = bmax;

    c = &common;
    CHOLMOD(start)(c);
    A = CHOLMOD(allocate_sparse)(M, N, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, c);
    c_float *Ax;
    c_int *Ai, *Ap;
    Ax = A->x;
    Ap = A->p;
    Ai = A->i;
    Ax[0] = 1.0; Ax[1] = 3.0; Ax[2] = 5.0; Ax[3] = 2.0; Ax[4] = 4.0;  
    Ap[0] = 0; Ap[1] = 3; Ap[2] = 5;
    Ai[0] = 0; Ai[1] = 1; Ai[2] = 2; Ai[3] = 0; Ai[4] = 1;

    Q = CHOLMOD(allocate_sparse)(N, N, QNZMAX, TRUE, TRUE, -1, CHOLMOD_REAL, c);
    c_float *Qx;
    c_int *Qi, *Qp;
    Qx = Q->x;
    Qp = Q->p;
    Qi = Q->i;
    Qx[0] = 1.0; Qx[1] = -1.0; Qx[2] = -1.0; Qx[3] = 2.0; 
    Qp[0] = 0; Qp[1] = 2; Qp[2] = 4; 
    Qi[0] = 0; Qi[1] = 1; Qi[2] = 0; Qi[3] = 1; 

    data->A = A;
    data->Q = Q;
    
    // Setup workspace
    work = qpalm_setup(data, settings);

    CHOLMOD(finish)(c);

    c_free(data);
    c_free(settings);
}

void cholmod_suite_teardown(void) {
    qpalm_cleanup(work);

    CHOLMOD(start)(c);
    CHOLMOD(free_sparse)(&Q, c);
    CHOLMOD(free_sparse)(&A, c);
    CHOLMOD(finish)(c);

}

void cholmod_test_setup(void) {
    work->Qd[0] = 1.1; work->Qd[1] = -0.5;
    work->Ad[0] = 1.1; work->Ad[1] = -0.5; work->Ad[2] = 20;
    CHOLMOD(start)(&common);
}

void cholmod_test_teardown(void) {
    CHOLMOD(finish)(&common);
}

MU_TEST(test_mat_vec){
    mat_vec(A, work->chol->Qd, work->chol->Ad, c);
    mu_assert_double_eq(work->Ad[0], 0.1, TOL);
    mu_assert_double_eq(work->Ad[1], 1.3, TOL);
    mu_assert_double_eq(work->Ad[2], 5.5, TOL);

    mat_vec(Q, work->chol->Qd, work->chol->Qd, c);
    mu_assert_double_eq(work->Qd[0], 1.6, TOL);
    mu_assert_double_eq(work->Qd[1], -2.1, TOL);
    work->Qd[0] = 1.1; work->Qd[1] = -0.5;

    mat_tpose_vec(Q, work->chol->Qd, work->chol->Qd, c);
    mu_assert_double_eq(work->Qd[0], 1.6, TOL);
    mu_assert_double_eq(work->Qd[1], -2.1, TOL);

}

MU_TEST(test_mat_tpose_vec){
    mat_tpose_vec(A, work->chol->Ad, work->chol->Qd, c);
    mu_assert_double_eq(work->Qd[0], 99.6, TOL);
    mu_assert_double_eq(work->Qd[1], 0.2, TOL);
}

MU_TEST(test_mat_inf_norm_cols){
    mat_inf_norm_cols(A, work->D_temp);
    mu_assert_double_eq(work->D_temp[0], 5.0, TOL);
    mu_assert_double_eq(work->D_temp[1], 4.0, TOL);
}

MU_TEST(test_mat_inf_norm_rows){
    mat_inf_norm_rows(A, work->E_temp);
    mu_assert_double_eq(work->E_temp[0], 2.0, TOL);
    mu_assert_double_eq(work->E_temp[1], 4.0, TOL);
    mu_assert_double_eq(work->E_temp[2], 5.0, TOL);
}
MU_TEST(test_ldlchol){
    // without proximal
    work->settings->proximal = FALSE;
    work->dphi[0] = -1.0; work->dphi[1] = -2.0; //this is -rhs
    ldlchol(Q, work, c);
    ldlsolveLD_neg_dphi(work, c);
    mu_assert_double_eq(work->d[0], 4.0, TOL);
    mu_assert_double_eq(work->d[1], 3.0, TOL);

    // with proximal
    work->settings->proximal = TRUE;
    work->gamma = 1e3;
    ldlchol(Q, work, c);
    ldlsolveLD_neg_dphi(work, c);
    mu_assert_double_eq(work->d[0], 3.989028924198480, TOL);
    mu_assert_double_eq(work->d[1], 2.993017953122679, TOL);
}

MU_TEST_SUITE(suite_cholmod) {
    MU_SUITE_CONFIGURE(cholmod_suite_setup, cholmod_suite_teardown, cholmod_test_setup, cholmod_test_teardown);

    MU_RUN_TEST(test_mat_vec);
    MU_RUN_TEST(test_mat_tpose_vec);
    MU_RUN_TEST(test_mat_inf_norm_cols);
    MU_RUN_TEST(test_mat_inf_norm_rows);
    MU_RUN_TEST(test_ldlchol);
    
}