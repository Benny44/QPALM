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
            if (sigma_temp <= 1e8) { //TODO make sigma_max a setting
                if (work->sigma[k] != sigma_temp) {
                    sigma_changed[work->nb_sigma_changed] = k;
                    work->nb_sigma_changed++;
                }               
                work->sigma[k] = sigma_temp;
                mult_factor = c_sqrt(mult_factor);
                work->sqrt_sigma[k] = mult_factor * work->sqrt_sigma[k];
                At_scalex[k] = mult_factor;
            } else {
                if (work->sigma[k] != 1e8) {
                    sigma_changed[work->nb_sigma_changed] = k;
                    work->nb_sigma_changed++;
                } 
                work->sigma[k] = 1e8;
                At_scalex[k] = 1e4 / work->sqrt_sigma[k];
                work->sqrt_sigma[k] = 1e4;
            }
        } else {
            At_scalex[k] = 1.0;
        }
    }

    CHOLMOD(scale)(work->chol->At_scale, CHOLMOD_COL, work->chol->At_sqrt_sigma, &work->chol->c);


    if ((work->settings->proximal && work->gamma != work->settings->gamma_max) || (work->nb_sigma_changed > MAX_RANK_UPDATE)) {
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