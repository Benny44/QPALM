/**
 * @file nonconvex.c
 * @author Ben Hermans
 * @brief Routines to deal with nonconvex QPs.
 * @details The functions in this file serve to set up QPALM for a nonconvex QP. The main routine in 
 * this file computes the minimum eigenvalue of a square matrix, based on lobpcg \cite knyazev2001toward. Furthermore, 
 * some setting updates are performed. In addition, the spectrum of a matrix can be upper bounded using 
 * Gershgorin's circle theorem, which is used in the gamma_boost routine in iteration.c.
 */

#include "nonconvex.h"
#include "types.h"
#include "constants.h"
#include "global_opts.h"
#include "lin_alg.h"
#include "util.h"
#ifdef MATLAB /* For the mexinterface, call lapack included in matlab */
    #include "lapack.h"
#else
    #include <lapacke.h>
#endif

/*TODO: make this a setting */
#define LOBPCG_TOL 1e-5 /**< Tolerance on the infinity norm of the residual in lobpcg. */ 




static c_float lobpcg(QPALMWorkspace *work, c_float *x, solver_common *c) {
    c_float lambda, norm_w;
    size_t i;

    size_t n = work->data->n;
    solver_sparse* A = work->data->Q;
    // size_t m = work->data->m;

    /*Current guess of the eigenvector */
    if (x == NULL) {
        x = work->d; 
        /* Initialize eigenvector randomly. */
        for (i = 0; i < n; i++) {
            x[i] = (c_float) rand()/RAND_MAX;
        }
        vec_self_mult_scalar(x, 1.0/vec_norm_two(x, n), n);
    } else {
        /* NB: Assume x is already normalized */
        prea_vec_copy(x, work->d, n);
        x = work->d;
    }
    solver_dense *x_chol = work->solver->d;

    
    c_float *Ax = work->Qd;
    solver_dense * Ax_chol = work->solver->Qd;
    mat_vec(A, x_chol , Ax_chol, c);
    lambda = vec_prod(x, Ax, n);

    /*Current residual, Ax - lambda*x */
    c_float *w = work->neg_dphi; 
    solver_dense * w_chol = work->solver->neg_dphi;
    c_float *Aw = work->Atyh;
    solver_dense * Aw_chol = work->solver->Atyh;

    /* Conjugate gradient direction */
    c_float *p = work->temp_n; 
    c_float *Ap = work->xx0;
    c_float p_norm_inv;

    /* Compressed system B = [x, w, p]'*Q*[x, w, p] */
    c_float B[3][3]; 
    c_float C[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; /* C = [1 0 xp; 0 1 wp; xp wp 1] Takes into account that p is not orthonormal with x, w */
    c_float lambda_B[3]; /* The eigenvalues of B */
    c_float *y; /* Eigenvector corresponding to min(lambda_B) */
    c_float xAw, wAw, xAp, wAp, pAp, xp, wp;

    /* Compute residual and make it orthonormal wrt x */
    vec_add_scaled(Ax, x, w, -lambda, n);
    vec_add_scaled(w, x, w, -vec_prod(x, w, n), n);
    vec_self_mult_scalar(w, 1.0/vec_norm_two(w, n), n);
    mat_vec(A, w_chol, Aw_chol, c);
    xAw = vec_prod(Aw, x, n);
    wAw = vec_prod(Aw, w, n);

    /* In the first compressed system, there is no p yet, so it is 2 by 2 */
    c_float B_init[2][2] = {{lambda, xAw}, {xAw, wAw}};
    c_float lambda_init[2];

    /* Lapack variables */
    char jobz = 'V';
    char uplo = 'L';

    /* Solve eigenvalue problem */
    #ifdef MATLAB
        double lapack_work[10];
        long int info = 0, dim = 2, lwork = 10, itype = 1;
        dsyev(&jobz, &uplo, &dim, *B_init, &dim, lambda_init, lapack_work, &lwork, &info);
    #else
        int info = 0, dim = 2, itype = 1;
        info = LAPACKE_dsyev(LAPACK_COL_MAJOR, jobz, uplo, dim, *B_init, dim, lambda_init);
    #endif
    lambda = lambda_init[0];
    y = B_init[0];

    /* Compute first p */
    vec_mult_scalar(w, y[1], p, n);
    vec_mult_scalar(Aw, y[1], Ap, n);
    vec_add_scaled(p, x, x, y[0], n);
    vec_add_scaled(Ap, Ax, Ax, y[0], n);
    
    dim = 3; /* From now on, the dimension of the eigenproblem to solve will be 3 */
    size_t max_iter = 1000; /*TODO: make this a setting */
    for (i = 0; i < max_iter; i++) {

        /* Update w */
        vec_add_scaled(Ax, x, w, -lambda, n);
        /* Note: check inf norm because it is cheaper */
        if (vec_norm_inf(w, n) < LOBPCG_TOL) {
            norm_w = vec_norm_two(w, n);
            lambda -= c_sqrt(2)*norm_w + 1e-6; /* Theoretical bound on the eigenvalue */
            if (n <= 3) lambda -= 1e-6; /* If n <= 3, we should have the exact eigenvalue, hence we subtract a small value */
            break;
        } 
        vec_add_scaled(w, x, w, -vec_prod(x, w, n), n);
        vec_self_mult_scalar(w, 1.0/vec_norm_two(w, n), n);
        mat_vec(A, w_chol, Aw_chol, c);
        xAw = vec_prod(Ax, w, n);
        wAw = vec_prod(w, Aw, n);

        /* Normalize p */
        p_norm_inv = 1.0/vec_norm_two(p, n);
        vec_self_mult_scalar(p, p_norm_inv, n);
        vec_self_mult_scalar(Ap, p_norm_inv, n);

        /* Compress the system */
        xAp = vec_prod(Ax, p, n);
        wAp = vec_prod(Aw, p, n);
        pAp = vec_prod(Ap, p, n);
        xp = vec_prod(x, p, n);
        wp = vec_prod(w, p, n);

        B[0][0] = lambda; B[0][1] = xAw; B[0][2] = xAp; 
        B[1][0] = xAw;    B[1][1] = wAw; B[1][2] = wAp; 
        B[2][0] = xAp;    B[2][1] = wAp; B[2][2] = pAp;

        C[0][2] = xp; C[1][2] = wp; C[2][0] = xp; C[2][1] = wp; 
        C[2][2] = 1.0; /* The dsygv routine might override this element, therefore we reset it here.*/

        /* Solve eigenproblem B*x = lambda*C*x */
        #ifdef MATLAB
            dsygv(&itype, &jobz, &uplo, &dim, *B, &dim, *C, &dim, lambda_B, lapack_work, &lwork, &info);
        #else
            info = LAPACKE_dsygv(LAPACK_COL_MAJOR, itype, 'V', 'L', dim, *B, dim, *C, dim, lambda_B);
        #endif
        lambda = lambda_B[0];
        y = B[0];

        /* Update p and x */
        vec_mult_add_scaled(p, w, y[2], y[1], n);
        vec_mult_add_scaled(Ap, Aw, y[2], y[1], n);
        vec_mult_add_scaled(x, p, y[0], 1, n);
        vec_mult_add_scaled(Ax, Ap, y[0], 1, n);

    }
    
    //TODO: Implement error handling here
    return lambda;

}


void set_settings_nonconvex(QPALMWorkspace *work, solver_common *c){
    c_float lambda;
    lambda = lobpcg(work, NULL, c);
    if (lambda < 0) {
        work->settings->proximal = TRUE;
        work->settings->gamma_init = 1/c_absval(lambda);
        work->settings->gamma_max = work->settings->gamma_init;
        work->gamma_maxed = TRUE;
    } else
    {
        work->settings->nonconvex = FALSE;
    }
}

c_float gershgorin_max(solver_sparse* M, c_float *center, c_float *radius){
    /* NB: Assume M is symmetric, so Gershgorin may be performed along the columns as well. */
    c_float ub_eig;
    c_float *Mx = M->x; c_int *Mi = M->i; c_int *Mp = M->p;
    c_int row, i, j, ncol = (c_int)M->ncol;
    
    for (i=0; i < ncol; i++) {
        center[i] = 0.0;
        radius[i] = 0.0;
        for (j = Mp[i]; j < Mp[i+1]; j++) {
            row = Mi[j];
            if (row == i) {
                center[i] = Mx[j];
            } else {
                radius[i] += c_absval(Mx[j]);
            }    
        }
        if (i==0) {
            ub_eig = center[i] + radius[i];
        } else {
            ub_eig = c_max(ub_eig, (center[i] + radius[i]));
        }
    }

    return ub_eig;
}

