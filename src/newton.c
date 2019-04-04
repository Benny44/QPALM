#include "newton.h"
#include "lin_alg.h"
#include "cholmod.h"
#include <stdio.h>

void newton_set_direction(QPALMWorkspace *work) {

    set_active_constraints(work);
    work->chol->nb_active_constraints = 0;
    if (work->chol->nb_active_constraints) {

    } else {
        //LD = ldlchol(Q+1/gamma*I)
        double beta [2] = {1.0/work->settings->gamma,0};
        work->chol->LD = cholmod_analyze (work->data->Q, &work->chol->c) ;
        cholmod_factorize_p (work->data->Q, beta, NULL, 0, work->chol->LD, &work->chol->c);
        //d = ldlsolve(LD, -dphi)
        prea_vec_copy(work->dphi, work->neg_dphi, work->data->n);
        vec_mult_scalar(work->neg_dphi, -1, work->data->n);
        work->chol->d = cholmod_solve (CHOLMOD_LDLt, work->chol->LD, work->chol->neg_dphi, &work->chol->c);
        work->d = work->chol->d->x;
    }

}

void set_active_constraints(QPALMWorkspace *work) {
    work->chol->nb_active_constraints = 0;
    for (int i = 0; i < work->data->m; i++) {
        if ((work->Axys[i] < work->data->bmin[i]) || ((work->Axys[i] > work->data->bmax[i]))){
            work->chol->active_constraints[i] = 1;
            work->chol->nb_active_constraints++;
        } else {
            work->chol->active_constraints[i] = 0;
        }         
    }
}