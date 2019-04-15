#include "cholmod_interface.h"
#include "global_opts.h"
#include "qpalm.h"
#include "constants.h"
#include <CUnit/CUnit.h>

#define N 2
#define M 3
#define ANZMAX 5
#define QNZMAX 4
#define TOL 1e-8

QPALMWorkspace *work; // Workspace
cholmod_sparse *A; // MxN matrix
cholmod_sparse *Q; //NxN symmetric matrix

void cholmod_qp_setup(void) {
    QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;

    QPALMData *data    = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;
    c_float q[N] = {1, -2};
    data->q = q;
    c_float bmin[M] = {-5, -10, -2};
    c_float bmax[M] = {5, 10, 3};
    data->bmin = bmin;
    data->bmax = bmax;

    cholmod_common c;
    CHOLMOD(start)(&c);
    A = CHOLMOD(allocate_sparse)(M, N, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    c_float *Ax;
    c_int *Ai, *Ap;
    Ax = A->x;
    Ap = A->p;
    Ai = A->i;
    Ax[0] = 1.0; Ax[1] = 3.0; Ax[2] = 5.0; Ax[3] = 2.0; Ax[4] = 4.0;  
    Ap[0] = 0; Ap[1] = 3; Ap[2] = 5;
    Ai[0] = 0; Ai[1] = 1; Ai[2] = 2; Ai[3] = 0; Ai[4] = 1;

    Q = CHOLMOD(allocate_sparse)(N, N, QNZMAX, TRUE, TRUE, -1, CHOLMOD_REAL, &c);
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
    CHOLMOD(finish)(&c);
    // Setup workspace
    work = qpalm_setup(data, settings, &c);

    c_free(data);
    c_free(settings);
}

void cholmod_qp_teardown(void) {
    qpalm_cleanup(work);

    cholmod_common c;
    CHOLMOD(start)(&c);
    CHOLMOD(free_sparse)(&Q, &c);
    CHOLMOD(free_sparse)(&A, &c);
    CHOLMOD(finish)(&c);
}

void cholmod_set_QdAd(void) {
    work->Qd[0] = 1.1; work->Qd[1] = -0.5;
    work->Ad[0] = 1.1; work->Ad[1] = -0.5; work->Ad[2] = 20;
}

void test_mat_vec(void){
    mat_vec(A, work->chol->Qd, work->chol->Ad, &work->chol->c);
    CU_ASSERT_DOUBLE_EQUAL(work->Ad[0], 0.1, TOL);
    CU_ASSERT_DOUBLE_EQUAL(work->Ad[1], 1.3, TOL);
    CU_ASSERT_DOUBLE_EQUAL(work->Ad[2], 5.5, TOL);

    mat_vec(Q, work->chol->Qd, work->chol->Qd, &work->chol->c);
    CU_ASSERT_DOUBLE_EQUAL(work->Qd[0], 1.6, TOL);
    CU_ASSERT_DOUBLE_EQUAL(work->Qd[1], -2.1, TOL);
}

void test_mat_tpose_vec(void){
    mat_tpose_vec(A, work->chol->Ad, work->chol->Qd, &work->chol->c);
    CU_ASSERT_DOUBLE_EQUAL(work->Qd[0], 99.6, TOL);
    CU_ASSERT_DOUBLE_EQUAL(work->Qd[1], 0.2, TOL);
}

void test_mat_inf_norm_cols(void){
    mat_inf_norm_cols(A, work->D_temp);
    CU_ASSERT_DOUBLE_EQUAL(work->D_temp[0], 5.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(work->D_temp[1], 4.0, TOL);
}

void test_mat_inf_norm_rows(void){
    mat_inf_norm_rows(A, work->E_temp);
    CU_ASSERT_DOUBLE_EQUAL(work->E_temp[0], 2.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(work->E_temp[1], 4.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(work->E_temp[2], 5.0, TOL);
}
void test_ldlchol(void){
    // without proximal
    work->settings->proximal = FALSE;
    work->dphi[0] = -1.0; work->dphi[1] = -2.0; //this is -rhs
    ldlchol(Q, work);
    ldlsolveLD_neg_dphi(work);
    CU_ASSERT_DOUBLE_EQUAL(work->d[0], 4.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(work->d[1], 3.0, TOL);

    // with proximal
    work->settings->proximal = TRUE;
    work->settings->gamma = 1e3;
    ldlchol(Q, work);
    ldlsolveLD_neg_dphi(work);
    CU_ASSERT_DOUBLE_EQUAL(work->d[0], 3.989028924198480, TOL);
    CU_ASSERT_DOUBLE_EQUAL(work->d[1], 2.993017953122679, TOL);
}