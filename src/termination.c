#include "termination.h"
#include "lin_alg.h"
#include "constants.h"
#include "util.h"

#include <stdio.h>

c_int check_termination(QPALMWorkspace *work) {
    calculate_residuals_and_tolerances(work);

    if (is_solved(work)) {
        update_status(work->info, QPALM_SOLVED);
        store_solution(work);
        return 1;
    } else if (is_primal_infeasible(work)) {
        update_status(work->info, QPALM_PRIMAL_INFEASIBLE);
        if (work->settings->scaling) {
            vec_ew_prod(work->scaling->E, work->delta_y, work->delta_y, work->data->m);
            vec_mult_scalar(work->y, work->scaling->cinv, work->data->m);
        }
        return 1;
    } else if (is_dual_infeasible(work)) {
        update_status(work->info, QPALM_DUAL_INFEASIBLE);
        if (work->settings->scaling) {
            vec_ew_prod(work->scaling->D, work->delta_x, work->delta_x, work->data->n);
        }
        return 1;
    } else {
        return 0;
    }   
}

void calculate_residuals_and_tolerances(QPALMWorkspace *work) {
    calculate_primal_residual(work);
    calculate_dual_residual(work);
    calculate_primal_tolerance(work);
    calculate_dual_tolerances(work);
}

void calculate_primal_residual(QPALMWorkspace *work) {
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Einv, work->pri_res, work->temp_m, work->data->m);
        work->info->pri_res_norm = vec_norm_inf(work->temp_m, work->data->m);
    } else {
        work->info->pri_res_norm = vec_norm_inf(work->pri_res, work->data->m);
    }
}

void calculate_dual_residual(QPALMWorkspace *work) {
    if (work->settings->scaling) {
        if (work->settings->proximal) {
            vec_add_scaled(work->dphi, work->xx0, work->temp_n, -1/work->settings->gamma, work->data->n);
            vec_ew_prod(work->scaling->Dinv, work->temp_n, work->temp_n, work->data->n);
            work->info->dua_res_norm = vec_norm_inf(work->temp_n, work->data->n);
            vec_ew_prod(work->scaling->Dinv, work->dphi, work->temp_n, work->data->n);
            work->info->dua2_res_norm = vec_norm_inf(work->temp_n, work->data->n);
        } else {
            vec_ew_prod(work->scaling->Dinv, work->dphi, work->temp_n, work->data->n);
            work->info->dua_res_norm = vec_norm_inf(work->temp_n, work->data->n);
            work->info->dua2_res_norm = work->info->dua_res_norm;
        } 
        work->info->dua_res_norm *= work->scaling->cinv;
        work->info->dua2_res_norm *= work->scaling->cinv;
            
    } else {
        if (work->settings->proximal) {
            vec_add_scaled(work->dphi, work->xx0, work->temp_n, -1/work->settings->gamma, work->data->n);
            work->info->dua_res_norm = vec_norm_inf(work->temp_n, work->data->n);
            work->info->dua2_res_norm = vec_norm_inf(work->dphi, work->data->n);
        } else {
            work->info->dua_res_norm = vec_norm_inf(work->dphi, work->data->n);
            work->info->dua2_res_norm = work->info->dua_res_norm;
        }
    }
}

void calculate_primal_tolerance(QPALMWorkspace *work) {
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Einv, work->Ax, work->temp_2m, work->data->m);
        vec_ew_prod(work->scaling->Einv, work->z, work->temp_2m + work->data->m, work->data->m);
        work->eps_pri =  work->settings->eps_abs + work->settings->eps_rel*vec_norm_inf(work->temp_2m, work->data->m);                  
    } else {
        work->eps_pri =  work->settings->eps_abs + work->settings->eps_rel*c_max(
                                vec_norm_inf(work->Ax, work->data->m),
                                vec_norm_inf(work->z, work->data->m));
    }
}

void calculate_dual_tolerances(QPALMWorkspace *work) {
    c_float norm_DinvQx, norm_Dinvq, norm_DinvAtyh, max_norm;
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Dinv, work->Qx, work->temp_n, work->data->n);
        norm_DinvQx = vec_norm_inf(work->temp_n, work->data->n);
        vec_ew_prod(work->scaling->Dinv, work->data->q, work->temp_n, work->data->n);
        norm_Dinvq = vec_norm_inf(work->temp_n, work->data->n);
        vec_ew_prod(work->scaling->Dinv, work->Atyh, work->temp_n, work->data->n);
        norm_DinvAtyh = vec_norm_inf(work->temp_n, work->data->n);
    } else {
        norm_DinvQx = vec_norm_inf(work->Qx, work->data->n);
        norm_Dinvq = vec_norm_inf(work->data->q, work->data->n);
        norm_DinvAtyh = vec_norm_inf(work->Atyh, work->data->n);
    }
    
    max_norm = c_max(norm_DinvQx, c_max(norm_Dinvq, norm_DinvAtyh));
    if (work->settings->scaling) max_norm *= work->scaling->cinv;

    work->eps_dua = work->settings->eps_abs + work->settings->eps_rel*max_norm;
    work->eps_dua_in = work->settings->eps_abs_in + work->settings->eps_rel_in*max_norm;
}

c_int is_solved(QPALMWorkspace *work) {
    return (work->info->pri_res_norm < work->eps_pri) && 
            (work->info->dua_res_norm < work->eps_dua);
}

c_int is_primal_infeasible(QPALMWorkspace *work) {
    c_float eps_pinf_norm_Edy;

    //dy = yh-y
    vec_add_scaled(work->yh, work->y, work->delta_y, -1, work->data->m);
    if (work->settings->scaling) {
        //Edy = E.*dy
        vec_ew_prod(work->scaling->E, work->delta_y, work->temp_m, work->data->m);
        eps_pinf_norm_Edy = work->settings->eps_prim_inf*
        vec_norm_inf(work->temp_m, work->data->m);
    } else {
        eps_pinf_norm_Edy = work->settings->eps_prim_inf*
        vec_norm_inf(work->delta_y, work->data->m);
    }
     
    if (eps_pinf_norm_Edy == 0) { //dy == 0
        return 0;
    }

    //mat_tpose_vec(work->data->A, work->delta_y, work->Atdelta_y, 0, 0);
    vec_add_scaled(work->Atyh, work->Aty, work->Atdelta_y, -1, work->data->n);
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Dinv, work->Atdelta_y, work->Atdelta_y, work->data->n);
    }

    //out_of_bounds = bmax'*max(dy,0) + bmin'*min(dy,0)
    c_float out_of_bounds = 0;
    for(c_int i=0; i < work->data->m; i++) {
        out_of_bounds += work->data->bmax[i]*c_max(work->delta_y[i], 0)
            + work->data->bmin[i]*c_min(work->delta_y[i], 0);
    }

    return (vec_norm_inf(work->Atdelta_y, work->data->n) <= eps_pinf_norm_Edy)
        && (out_of_bounds <= -eps_pinf_norm_Edy);
    
}

c_int is_dual_infeasible(QPALMWorkspace *work) {
    c_float eps_dinf_norm_Ddx;

    //dx = x-x_prev
    vec_add_scaled(work->x, work->x_prev, work->delta_x, -1, work->data->n);
    if (work->settings->scaling) {
        //D*dx
        vec_ew_prod(work->scaling->D, work->delta_x, work->temp_n, work->data->n);
        eps_dinf_norm_Ddx = work->settings->eps_dual_inf*
                                vec_norm_inf(work->temp_n, work->data->n);
    } else {
        eps_dinf_norm_Ddx = work->settings->eps_dual_inf*
                                vec_norm_inf(work->delta_x, work->data->n);
    }
    
    if (eps_dinf_norm_Ddx == 0) { //dx == 0
        return 0;
    }

    //NB Adx = work->Ad (= tau*Ad of the previous iteration)
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Einv, work->Ad, work->Adelta_x, work->data->m);
        for (c_int k = 0; k < work->data->m; k++) {
            if ((work->data->bmax[k] < QPALM_INFTY && work->Adelta_x[k] >= eps_dinf_norm_Ddx)
                || (work->data->bmin[k] > -QPALM_INFTY && work->Adelta_x[k] <= -eps_dinf_norm_Ddx)) {
                return 0;
            }
        }      
    } else {
        for (c_int k = 0; k < work->data->m; k++) {
            if ((work->data->bmax[k] < QPALM_INFTY && work->Ad[k] >= eps_dinf_norm_Ddx)
                || (work->data->bmin[k] > -QPALM_INFTY && work->Ad[k] <= -eps_dinf_norm_Ddx)) {
                return 0;
            }
        }
    }
    //NB Qdx = work->Qd (= tau*Qd of the previous iteration)
    if (work->settings->scaling) {
        vec_ew_prod(work->scaling->Dinv, work->Qd, work->temp_n, work->data->n);
        return (vec_norm_inf(work->temp_n, work->data->n) <= work->scaling->c*eps_dinf_norm_Ddx)
            && (vec_prod(work->data->q, work->delta_x, work->data->n) <= -work->scaling->c*eps_dinf_norm_Ddx);
    } else {
        return (vec_norm_inf(work->Qd, work->data->n) <= eps_dinf_norm_Ddx)
            && (vec_prod(work->data->q, work->delta_x, work->data->n) <= -eps_dinf_norm_Ddx);
    }
    
}

void store_solution(QPALMWorkspace *work) {
    if (work->settings->scaling) {
        vec_ew_prod(work->x, work->scaling->D, work->solution->x, work->data->n);
        vec_ew_prod(work->yh, work->scaling->E, work->solution->y, work->data->m);
        vec_mult_scalar(work->yh, work->scaling->cinv, work->data->m);
    } else {
        prea_vec_copy(work->x, work->solution->x, work->data->n);
        prea_vec_copy(work->yh, work->solution->y, work->data->m);
    }
    
}

c_int check_subproblem_termination(QPALMWorkspace *work) {
    return (work->info->dua2_res_norm <= work->eps_dua_in);
}