#include "linesearch.h"
#include "lin_alg.h"
#include <stdlib.h> //for sorting

c_float exact_linesearch_newton(QPALMWorkspace *work){
    //Qd
    mat_vec_triu(work->data->Q, work->d, work->Qd);
    if (work->settings->proximal) {
        vec_add_scaled(work->Qd, work->d, work->Qd, 1/work->settings->gamma, work->data->n);
    }
    //Ad
    mat_vec(work->data->A, work->d, work->Ad);
    //eta = d'*Qd
    work->eta = vec_prod(work->d, work->Qd, work->data->n);
    //beta = d'*df
    work->beta = vec_prod(work->d, work->df, work->data->n);
    //sigma_sqrt = sqrt(sigma)
    vec_ew_sqrt(work->sigma, work->sigma_sqrt, work->data->m);

    //delta = [-sqrt(sigma).*Ad; sqrt(sigma).*Ad]
    vec_ew_prod(work->sigma_sqrt, work->Ad, work->temp_m, work->data->m);
    prea_vec_copy(work->temp_m, work->delta + work->data->m, work->data->m); //shifted copy
    vec_mult_scalar(work->temp_m, -1, work->data->m);
    prea_vec_copy(work->temp_m, work->delta, work->data->m); 
    //alpha = [(y+sigma.*(Ax-bmin))./sigma_sqrt; (-y+sigma.*(bmax-Ax))./sigma_sqrt]
    vec_add_scaled(work->Ax, work->data->bmin, work->temp_m, -1, work->data->m);
    vec_ew_prod(work->sigma, work->temp_m, work->temp_m, work->data->m);
    vec_add_scaled(work->y, work->temp_m, work->temp_m, 1, work->data->m);
    vec_ew_div(work->temp_m, work->sigma_sqrt, work->temp_m, work->data->m);
    prea_vec_copy(work->temp_m, work->alpha, work->data->m);
    vec_add_scaled(work->data->bmax, work->Ax, work->temp_m, -1, work->data->m);
    vec_ew_prod(work->sigma, work->temp_m, work->temp_m, work->data->m);
    vec_add_scaled(work->temp_m, work->y, work->temp_m, -1, work->data->m);
    vec_ew_div(work->temp_m, work->sigma_sqrt, work->temp_m, work->data->m);
    prea_vec_copy(work->temp_m, work->alpha + work->data->m, work->data->m); //shifted copy
    // s = alpha./delta
    vec_ew_div(work->alpha, work->delta, work->temp_2m, work->data->m*2);

    // index_P = delta > 0
    // index_J = delta < 0
    for (c_int i=0; i<work->data->m*2; i++){
        if (work->delta[i] == 0) {
            work->index_P[i] = 0;
            work->index_J[i] = 0;
        }
        else if (work->delta[i] > 0) {
            work->index_P[i] = 1;
            work->index_J[i] = 0;
        } else if (work->delta[i] < 0){
            work->index_P[i] = 0;
            work->index_J[i] = 1;
        }     
    }

    c_float intercept, slope, tau;
    tau = -QPALM_INFTY;
    //tau = max(s)
    for (c_int i = 0; i < work->data->m*2; i++) {
        if (work->temp_2m[i] > tau) {
            tau = work->temp_2m[i];
        }
    }

    //Precompute delta.^2 and delta.*alpha
    vec_ew_prod(work->delta, work->delta, work->delta2, work->data->m*2);
    vec_ew_prod(work->delta, work->alpha, work->delta_alpha, work->data->m*2);

    //Newton iterations
    int A = 0; 
    for (c_int k = 0; k < 100; k++) {
        intercept = work->beta;
        slope = work->eta;
        for (c_int i = 0; i < work->data->m*2; i++) {
            // printf(" %d\n", ((work->index_P[i] && (tau > work->temp_2m[i])) || (work->index_J[i] && (tau < work->temp_2m[i]))));
            // if ((work->index_P[i] && (tau > work->temp_2m[i])) || (work->index_J[i] && (tau < work->temp_2m[i]))) {
            //     intercept -= work->delta_alpha[i];
            //     slope     += work->delta2[i];
            // } 
            A = ((work->index_P[i] && (tau > work->temp_2m[i])) || (work->index_J[i] && (tau < work->temp_2m[i])));
            intercept -= A*work->delta_alpha[i];
            slope     += A*work->delta2[i];
        }
        if (c_absval(slope*tau+intercept) < 1e-12) {
            break;
        }
        tau = -intercept/slope;
    }
    return tau;
}

c_float armijo_linesearch(QPALMWorkspace *work) {
    //Qd
    mat_vec_triu(work->data->Q, work->d, work->Qd);
    if (work->settings->proximal) {
        vec_add_scaled(work->Qd, work->d, work->Qd, 1/work->settings->gamma, work->data->n);
    }
    //Ad
    mat_vec(work->data->A, work->d, work->Ad);

    c_float c1 = 1e-4; //TODO: make this a setting
    c_float rhs, lhs, dQd, ddf, dist2;
    dQd = vec_prod(work->d, work->Qd, work->data->n);
    ddf = vec_prod(work->d, work->df, work->data->n);
    rhs = c1*(ddf + vec_prod(work->Ad, work->yh, work->data->m));
    
    c_int i;
    dist2 = 0;
    for (i = 0; i < work->data->m; i++) {
        dist2 += work->sigma[i]*(work->Axys[i] - work->z[i])*(work->Axys[i] - work->z[i]);
    }

    c_float tau = 2.0*work->settings->tau_init;
    do {
        tau *= 0.5;
        lhs = 0.5*tau*tau*dQd + tau*ddf - 0.5*dist2;
        for (i = 0; i < work->data->m; i++) {
            if (work->Axys[i] + tau*work->Ad[i] > work->data->bmax[i]) {
                lhs += 0.5*work->sigma[i]*(work->Axys[i] + tau*work->Ad[i] - work->data->bmax[i])
                        *(work->Axys[i] + tau*work->Ad[i] - work->data->bmax[i]);
            } else if (work->Axys[i] + tau*work->Ad[i] < work->data->bmin[i]) {
                lhs += 0.5*work->sigma[i]*(work->Axys[i] + tau*work->Ad[i] - work->data->bmin[i])
                        *(work->Axys[i] + tau*work->Ad[i] - work->data->bmin[i]);
            }
        }
        // printf("lhs: %f\n", lhs);
    } while (lhs > tau*rhs);

    // printf("lhs: %f, rhs: %f\n", lhs, rhs);
    // printf("tau: %f\n", tau);
    return tau;    
}

c_float exact_linesearch(QPALMWorkspace *work) {

    //Qd
    mat_vec_triu(work->data->Q, work->d, work->Qd);
    if (work->settings->proximal) {
        vec_add_scaled(work->Qd, work->d, work->Qd, 1/work->settings->gamma, work->data->n);
    }
    //Ad
    mat_vec(work->data->A, work->d, work->Ad);
    //eta = d'*Qd
    work->eta = vec_prod(work->d, work->Qd, work->data->n);
    //beta = d'*df
    work->beta = vec_prod(work->d, work->df, work->data->n);
    //sigma_sqrt = sqrt(sigma)
    vec_ew_sqrt(work->sigma, work->sigma_sqrt, work->data->m);

    //delta = [-sqrt(sigma).*Ad; sqrt(sigma).*Ad]
    vec_ew_prod(work->sigma_sqrt, work->Ad, work->temp_m, work->data->m);
    prea_vec_copy(work->temp_m, work->delta + work->data->m, work->data->m); //shifted copy
    vec_mult_scalar(work->temp_m, -1, work->data->m);
    prea_vec_copy(work->temp_m, work->delta, work->data->m); 
    //alpha = [(y+sigma.*(Ax-bmin))./sigma_sqrt; (-y+sigma.*(bmax-Ax))./sigma_sqrt]
    vec_add_scaled(work->Ax, work->data->bmin, work->temp_m, -1, work->data->m);
    vec_ew_prod(work->sigma, work->temp_m, work->temp_m, work->data->m);
    vec_add_scaled(work->y, work->temp_m, work->temp_m, 1, work->data->m);
    vec_ew_div(work->temp_m, work->sigma_sqrt, work->temp_m, work->data->m);
    prea_vec_copy(work->temp_m, work->alpha, work->data->m);
    vec_add_scaled(work->data->bmax, work->Ax, work->temp_m, -1, work->data->m);
    vec_ew_prod(work->sigma, work->temp_m, work->temp_m, work->data->m);
    vec_add_scaled(work->temp_m, work->y, work->temp_m, -1, work->data->m);
    vec_ew_div(work->temp_m, work->sigma_sqrt, work->temp_m, work->data->m);
    prea_vec_copy(work->temp_m, work->alpha + work->data->m, work->data->m); //shifted copy
    // s = alpha./delta
    vec_ew_div(work->alpha, work->delta, work->temp_2m, work->data->m*2);
    vec_array_copy(work->temp_2m, work->s, work->data->m*2);
    // index_L = s>0
    c_int nL = 0;
    for (c_int i=0; i<work->data->m*2; i++){
        if (work->temp_2m[i] > 0) {
            work->index_L[i] = 1;
            nL++;
        } else {
            work->index_L[i] = 0;
        }     
    };
    // printf("nL: %d\n", (int) nL);
    //s = s(indL)
    select_subsequence(work->s, work->s, work->index_L, work->data->m*2);

    // index_P = delta > 0
    for (c_int i=0; i<work->data->m*2; i++){
        if (work->delta[i] > 0) {
            work->index_P[i] = 1;
        } else {
            work->index_P[i] = 0;
        }     
    };
    // index_J = (P&~L)|(~P&L);
    for (c_int i=0; i<work->data->m*2; i++){
        if ((work->index_P[i] + work->index_L[i]) == 1) {
            work->index_J[i] = 1;
        } else {
            work->index_J[i] = 0;
        }     
    };

    // a = eta+delta(J)'*delta(J);
    // b = beta-delta(J)'*alpha(J);
    c_float a, b;
    a = work->eta + vec_prod_ind(work->delta, work->delta, work->index_J, work->data->m*2);
    b = work->beta - vec_prod_ind(work->delta, work->alpha, work->index_J, work->data->m*2);
    // return 0;
    //s = sort(s)
    // qsort(work->s, work->data->m*2, sizeof(array_element), compare);
    qsort(work->s, nL, sizeof(array_element), compare);
    
    if (nL == 0 || a*work->s[0].x+b > 0) {
        return -b/a; 
    }

    c_int i = 0;
    c_int iz;
    while (i < nL-1) {
        iz = work->s[i].i;
        if (work->index_P[iz]) {
            a = a + work->delta[iz]*work->delta[iz];
            b = b - work->delta[iz]*work->alpha[iz];
        } else {
            a = a - work->delta[iz]*work->delta[iz];
            b = b + work->delta[iz]*work->alpha[iz];
        }
        i++;
        if (a*work->s[i].x+b > 0) {
            return -b/a;
        }
    }
    iz = work->s[i].i;
    if (work->index_P[iz]) {
        a = a + work->delta[iz]*work->delta[iz];
        b = b - work->delta[iz]*work->alpha[iz];
    }
    return -b/a;

}

void vec_array_copy(c_float *a, array_element* b, c_int n) {
    c_int i;
    array_element ae;

    for (i=0; i < n; i++) {
        ae.x = a[i];
        ae.i = i;
        b[i] = ae; 
    }

}

void select_subsequence(const array_element *a, array_element *b, const c_int *L, c_int n) {
    c_int i;
    c_int nb_elements = 0;
    for (i=0; i < n; i++) {
        if (L[i]) {
            b[nb_elements] = a[i];
            nb_elements++;
        }
    }
}

c_float vec_prod_ind(const c_float *a, const c_float *b, const c_int *L, c_int n) {
  c_float prod = 0.0;
  c_int   i; // Index

  for (i = 0; i < n; i++) {
      if (L[i]) {
          prod += a[i] * b[i];
      }   
  }

  return prod;
}

int compare (const void * a, const void * b)
{
    c_float f = ((struct array_element*)a)->x;
    c_float s = ((struct array_element*)b)->x;
    if (f > s) return  1;
    if (f < s) return -1;
    return 0;

}