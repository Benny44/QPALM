#include "iteration.h"
#include "lin_alg.h"
#include "cholmod_interface.h"
#include "newton.h"
#include "linesearch.h"

void compute_residuals(QPALMWorkspace *work) {

    //Axys = Ax + y./sigma
    vec_ew_div(work->y, work->sigma, work->temp_m, work->data->m);
    vec_add_scaled(work->Ax, work->temp_m, work->Axys, 1, work->data->m);
    //z = min(max(Axys,bmin),bmax)
    vec_ew_mid_vec(work->Axys, work->data->bmin, work->data->bmax, work->z, work->data->m);
    //pri_res = Ax-z
    vec_add_scaled(work->Ax, work->z, work->pri_res, -1, work->data->m);
    //yh = y + pri_res.*sigma
    vec_ew_prod(work->pri_res, work->sigma, work->temp_m, work->data->m);
    vec_add_scaled(work->y, work->temp_m, work->yh, 1, work->data->m);
    //df = Qx + q
    vec_add_scaled(work->Qx, work->data->q, work->df, 1, work->data->n);
    
    if (work->settings->proximal) {
      //df = Qx + q +1/gamma*(x-x0)
      // NB work->Qx contains Qx+1/gamma*x
      vec_add_scaled(work->df, work->x0, work->df, -1/work->gamma, work->data->n);
    }
    // Atyh = A'*yh
    mat_tpose_vec(work->data->A, work->chol->yh, work->chol->Atyh, &work->chol->c);
    //dphi = df+Atyh
    vec_add_scaled(work->df, work->Atyh, work->dphi, 1, work->data->n);
}

void initialize_sigma(QPALMWorkspace *work) {
    // Compute initial sigma
    size_t n = work->data->n;
    size_t m = work->data->m;
    c_float f = 0.5*vec_prod(work->x, work->Qx, n) + vec_prod(work->data->q, work->x, n);
    vec_ew_mid_vec(work->Ax, work->data->bmin, work->data->bmax, work->temp_m, m);
    vec_add_scaled(work->Ax, work->temp_m, work->temp_m, -1, m);
    c_float dist2 = vec_prod(work->temp_m, work->temp_m, m);
    vec_set_scalar(work->sigma, c_max(1e-4, c_min(2e1*c_max(1,c_absval(f))/c_max(1,0.5*dist2),1e4)), m);

    // Set fields related to sigma
    vec_ew_sqrt(work->sigma, work->sqrt_sigma, m);
    c_float *At_scalex = work->chol->At_scale->x;
    prea_vec_copy(work->sqrt_sigma, At_scalex, m);
    if (work->chol->At_sqrt_sigma) 
      CHOLMOD(free_sparse)(&work->chol->At_sqrt_sigma, &work->chol->c);
    work->chol->At_sqrt_sigma = CHOLMOD(transpose)(work->data->A, 1, &work->chol->c);
    CHOLMOD(scale)(work->chol->At_scale, CHOLMOD_COL, work->chol->At_sqrt_sigma, &work->chol->c);

    work->sqrt_sigma_max = c_sqrt(work->settings->sigma_max);
}

void update_sigma(QPALMWorkspace* work) {
    work->nb_sigma_changed = 0;
    c_float *At_scalex = work->chol->At_scale->x;
    c_float pri_res_unscaled_norm = vec_norm_inf(work->pri_res, work->data->m);
    c_float sigma_temp, mult_factor;
    c_int *sigma_changed = work->chol->enter;
    size_t k;
    for (k = 0; k < work->data->m; k++) {
        if ((c_absval(work->pri_res[k]) > work->settings->theta*c_absval(work->pri_res_in[k])) && work->chol->active_constraints[k]) {
            mult_factor = c_max(1.0, work->settings->delta * c_absval(work->pri_res[k]) / (pri_res_unscaled_norm + 1e-6));
            sigma_temp = mult_factor * work->sigma[k];
            if (sigma_temp <= work->settings->sigma_max) { 
                if (work->sigma[k] != sigma_temp) {
                    sigma_changed[work->nb_sigma_changed] = k;
                    work->nb_sigma_changed++;
                }               
                work->sigma[k] = sigma_temp;
                mult_factor = c_sqrt(mult_factor);
                work->sqrt_sigma[k] = mult_factor * work->sqrt_sigma[k];
                At_scalex[k] = mult_factor;
            } else {
                if (work->sigma[k] != work->settings->sigma_max) {
                    sigma_changed[work->nb_sigma_changed] = k;
                    work->nb_sigma_changed++;
                } 
                work->sigma[k] = work->settings->sigma_max;
                At_scalex[k] = work->sqrt_sigma_max / work->sqrt_sigma[k];
                work->sqrt_sigma[k] = work->sqrt_sigma_max;
            }
        } else {
            At_scalex[k] = 1.0;
        }
    }

    CHOLMOD(scale)(work->chol->At_scale, CHOLMOD_COL, work->chol->At_sqrt_sigma, &work->chol->c);
    


    if ((work->settings->proximal && work->gamma != work->settings->gamma_max) || (work->nb_sigma_changed > 0*MAX_RANK_UPDATE)) {
        work->chol->reset_newton = TRUE;
      } else if (work->nb_sigma_changed == 0){
        /* do nothing */
      } else {  
          ldlupdate_sigma_changed(work);
    }
}

void update_gamma(QPALMWorkspace *work) {
    c_float prev_gamma = work->gamma;
    work->gamma = c_min(work->gamma*work->settings->gamma_upd, work->settings->gamma_max);
    if (work->gamma != prev_gamma) work->chol->reset_newton = TRUE;
    prea_vec_copy(work->x, work->x0, work->data->n);
    vec_add_scaled(work->Qx, work->x, work->Qx, 1/work->gamma - 1/prev_gamma, work->data->n);
}

void update_primal_iterate(QPALMWorkspace *work) {
    newton_set_direction(work);

    work->tau = exact_linesearch(work);

    //x_prev = x
    prea_vec_copy(work->x, work->x_prev, work->data->n);
    //dphi_prev = dphi 
    prea_vec_copy(work->dphi, work->dphi_prev, work->data->n);
    //x = x+tau*d
    vec_add_scaled(work->x, work->d, work->x, work->tau, work->data->n);
    vec_mult_scalar(work->Qd, work->tau, work->data->n); //Qdx used in dua_infeas check
    vec_mult_scalar(work->Ad, work->tau, work->data->m); //Adx used in dua_infeas check
    vec_add_scaled(work->Qx, work->Qd, work->Qx, 1, work->data->n);
    vec_add_scaled(work->Ax, work->Ad, work->Ax, 1, work->data->m);
}

c_float compute_objective(QPALMWorkspace *work) {
    c_float objective = 0.0;
    size_t n = work->data->n;
    size_t i = 0; 

    if (work->settings->proximal) {
        if(n >= 4) {
            for (; i <= n-4; i+=4) {
                objective +=  (0.5*(work->Qx[i] - 1/work->gamma*work->x[i]) + work->data->q[i])*work->x[i] 
                            + (0.5*(work->Qx[i+1] - 1/work->gamma*work->x[i+1]) + work->data->q[i+1])*work->x[i+1]
                            + (0.5*(work->Qx[i+2] - 1/work->gamma*work->x[i+2]) + work->data->q[i+2])*work->x[i+2] 
                            + (0.5*(work->Qx[i+3] - 1/work->gamma*work->x[i+3]) + work->data->q[i+3])*work->x[i+3];
            }
        }
        for (; i < n; i++) {
            objective += (0.5*(work->Qx[i] - 1/work->gamma*work->x[i])+ work->data->q[i])*work->x[i];
        }
    } else {
        if(n >= 4) {
            for (; i <= n-4; i+=4) {
                objective += (0.5*work->Qx[i] + work->data->q[i])*work->x[i] 
                            + (0.5*work->Qx[i+1] + work->data->q[i+1])*work->x[i+1]
                            + (0.5*work->Qx[i+2] + work->data->q[i+2])*work->x[i+2] 
                            + (0.5*work->Qx[i+3] + work->data->q[i+3])*work->x[i+3];
            }
        }
        for (; i < n; i++) {
            objective += (0.5*work->Qx[i] + work->data->q[i])*work->x[i];
        }
    }

    if (work->settings->scaling) {
        objective *= work->scaling->cinv;
    }

    return objective;
}
