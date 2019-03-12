#include "lbfgs.h"
#include "lin_alg.h"
#include "global_opts.h"
#include <stdio.h>

void lbfgs_set_direction(QPALMWorkspace *work) {
    if (work->lbfgs->reset_lbfgs) {
        //d = -dphi
        prea_vec_copy(work->dphi, work->d, work->data->n);
        vec_mult_scalar(work->d, -1, work->data->n);
        work->lbfgs->reset_lbfgs = 0;
        work->lbfgs->curridx = 0;
        work->lbfgs->currmem = 0;
    } else {
        //s = x-x_prev
        vec_add_scaled(work->x, work->x_prev, work->lbfgs->s, -1, work->data->n);
        //y = dphi-dphi_prev
        vec_add_scaled(work->dphi, work->dphi_prev, work->lbfgs->y, -1, work->data->n);
        //ys = y'*s
        work->lbfgs->ys = vec_prod(work->lbfgs->y, work->lbfgs->s, work->data->n);

        if (work->lbfgs->ys >= c_max(1e-8*vec_norm_inf(work->lbfgs->s, work->data->n)*
            vec_norm_inf(work->lbfgs->y, work->data->n), 1e-12)) {
            
            lbfgs_update_buffers(work);
            //H = ys/(y'*y) 
            work->lbfgs->H = work->lbfgs->ys/vec_prod(work->lbfgs->y, work->lbfgs->y, work->data->n);
        }
        
        lbfgs_two_loop(work);
    }
}


void lbfgs_update_buffers(QPALMWorkspace *work) {
    c_int buffer_location;
    if (work->lbfgs->currmem == work->settings->memory) {
        buffer_location = mod(work->lbfgs->curridx, work->settings->memory)*work->data->n;
        prea_vec_copy(work->lbfgs->s, work->lbfgs->Sbuffer + buffer_location, work->data->n);
        prea_vec_copy(work->lbfgs->y, work->lbfgs->Ybuffer + buffer_location, work->data->n);
        work->lbfgs->YSbuffer[mod(work->lbfgs->curridx, work->settings->memory)] = work->lbfgs->ys;
        work->lbfgs->curridx = mod(work->lbfgs->curridx+1, work->settings->memory);
    } else {
        buffer_location = work->lbfgs->currmem*work->data->n;
        prea_vec_copy(work->lbfgs->s, work->lbfgs->Sbuffer + buffer_location, work->data->n);
        prea_vec_copy(work->lbfgs->y, work->lbfgs->Ybuffer + buffer_location, work->data->n);
        work->lbfgs->YSbuffer[work->lbfgs->currmem] = work->lbfgs->ys;
        work->lbfgs->currmem++;
        work->lbfgs->curridx++;
    }
}

void lbfgs_two_loop(QPALMWorkspace *work) {
    
    //alpha = zeros(1,mem)
    vec_set_scalar(work->lbfgs->alpha, 0, work->lbfgs->currmem);
    //q = -dphi
    prea_vec_copy(work->dphi, work->lbfgs->q, work->data->n);
    vec_mult_scalar(work->lbfgs->q, -1, work->data->n);

    c_int i = mod(work->lbfgs->curridx-1, work->settings->memory);
    c_int buffer_location = i*work->data->n;

    for (c_int k=0; k < work->lbfgs->currmem; k++) {
        work->lbfgs->alpha[i] = (1/work->lbfgs->YSbuffer[i])*
            vec_prod(work->lbfgs->Sbuffer + buffer_location, work->lbfgs->q, work->data->n);
        vec_add_scaled(work->lbfgs->q, work->lbfgs->Ybuffer + buffer_location, 
            work->lbfgs->q, -1*work->lbfgs->alpha[i], work->data->n);
        i = mod((i-1), work->settings->memory);
        buffer_location = i*work->data->n;
    }

    prea_vec_copy(work->lbfgs->q, work->d, work->data->n);
    vec_mult_scalar(work->d, work->lbfgs->H, work->data->n);

    i = mod(i+1, work->settings->memory);
    buffer_location = i*work->data->n;

    c_float beta;
    for (c_int k=0; k < work->lbfgs->currmem; k++) {
        beta = (1/work->lbfgs->YSbuffer[i])*
            vec_prod(work->lbfgs->Ybuffer + buffer_location, work->d, work->data->n);
        vec_add_scaled(work->d, work->lbfgs->Sbuffer + buffer_location, 
            work->d, work->lbfgs->alpha[i] - beta, work->data->n);
        i = mod(i+1, work->settings->memory);
        buffer_location = i*work->data->n;
    }

}