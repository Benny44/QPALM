/**
 * @file nonconvex.c
 * @author Ben Hermans
 * @brief Routines to deal with nonconvex QPs.
 * @details The functions in this file serve to set up QPALM for a nonconvex QP. The main routine in 
 * this file computes the minimum eigenvalue of a square matrix, based on power iterations. Furthermore, 
 * some setting updates are performed. 
 */

#include "nonconvex.h"
#include "types.h"
#include "constants.h"
#include "global_opts.h"
#include "lin_alg.h"

#define TOL 1e-3 /*TODO: make this a setting */

c_float minimum_eigenvalue_Q(QPALMWorkspace *work){
    c_float lambda;
    /* calculate largest (in absolute value) eigenvalue */
    lambda = power_iterations_Q(work, 0);
    if (lambda > 0) {
        /* calculate smallest eigenvalue by shifting */
        lambda += power_iterations_Q(work, -lambda); 
    }
    return lambda;
}

c_float power_iterations_Q(QPALMWorkspace *work, c_float gamma){
    /* NB: use neg_dphi and Qd as links to cholmod data structures */
    c_float lambda, lambda_prev, Qd_norm;

    size_t n = work->data->n;
    lambda = 0;
    lambda_prev = -QPALM_INFTY;
    vec_set_scalar(work->neg_dphi, 1.0, n);

    while(c_absval(lambda-lambda_prev) > TOL) {
        lambda_prev = lambda;
        /* Qd = (Q + gamma*I)*neg_dphi */ 
        mat_vec(work->data->Q, work->chol->neg_dphi, work->chol->Qd, &work->chol->c);
        vec_add_scaled(work->Qd, work->neg_dphi, work->Qd, gamma, n);

        Qd_norm = vec_norm_inf(work->Qd, n);
        if (Qd_norm == 0) { /*Q is zeros*/
            return TOL; /*zero after adjustments later*/
        }
        
        lambda = vec_prod(work->Qd, work->neg_dphi, n) / vec_prod(work->neg_dphi, work->neg_dphi, n);
        vec_mult_scalar(work->Qd, 1/Qd_norm, n);
        prea_vec_copy(work->Qd, work->neg_dphi, n);
    }

    return lambda;
}

void set_settings_nonconvex(QPALMWorkspace *work){
    c_float lambda;
    lambda = minimum_eigenvalue_Q(work);
    /*adjust for power iterations tolerance */
    lambda -= TOL; 
    if (lambda < 0) {
        work->settings->proximal = TRUE;
        work->settings->gamma_init = 1/c_absval(lambda);
        work->settings->gamma_max = work->settings->gamma_init;
    }
}


