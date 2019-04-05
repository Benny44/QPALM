#include "newton.h"
#include "lin_alg.h"
#include "cholmod.h"
#include <stdio.h>

void newton_set_direction(QPALMWorkspace *work) {

    set_active_constraints(work);
    if (work->chol->reset_newton && work->chol->nb_active_constraints) {
        work->chol->reset_newton = FALSE;
        //LD = ldlchol(Q+A(active_cnstrs,:)'*diag(sigma)*A(active_cnstrs,:));
        cholmod_sparse *AtsigmaA;
        cholmod_sparse *QAtsigmaA;
        size_t nb_active = 0;
        for (int i = 0; i < work->data->m; i++) {
            if (work->chol->active_constraints[i]){
                work->chol->enter[nb_active] = i;
                nb_active++;
            }
            
        }
        AtsigmaA = cholmod_aat(work->chol->At_sqrt_sigma, work->chol->enter, nb_active, TRUE, &work->chol->c);
        double one [2] = {1,0};
        QAtsigmaA = cholmod_add(work->data->Q, AtsigmaA, one, one, TRUE, FALSE, &work->chol->c);
        QAtsigmaA->stype = 1;
        c_float *Qx = QAtsigmaA->x;
        
        double beta [2] = {1.0/work->settings->gamma,0};
        work->chol->LD = cholmod_analyze (QAtsigmaA, &work->chol->c) ;
        cholmod_factorize_p (QAtsigmaA, beta, NULL, 0, work->chol->LD, &work->chol->c);

        cholmod_free_sparse(&AtsigmaA, &work->chol->c);
        cholmod_free_sparse(&QAtsigmaA, &work->chol->c);
        
    } else if (work->chol->nb_active_constraints) {
        set_entering_leaving_constraints(work);
        if(work->chol->nb_enter) {
            cholmod_sparse *Ae;
            Ae = cholmod_submatrix(work->chol->At_sqrt_sigma, NULL, -1, 
                                work->chol->enter, work->chol->nb_enter, TRUE, TRUE, &work->chol->c);
            //LD = ldlupdate(LD,Ae,'+');
            cholmod_updown(TRUE, Ae, work->chol->LD, &work->chol->c);
            cholmod_free_sparse(&Ae, &work->chol->c);
        }
        if(work->chol->nb_leave) {
            cholmod_sparse *Al;
            Al = cholmod_submatrix(work->chol->At_sqrt_sigma, NULL, -1, 
                                work->chol->leave, work->chol->nb_leave, TRUE, TRUE, &work->chol->c);
            //LD = ldlupdate(LD,Ae,'+');
            cholmod_updown(FALSE, Al, work->chol->LD, &work->chol->c);
            cholmod_free_sparse(&Al, &work->chol->c);
        }
    } else {
        //LD = ldlchol(Q+1/gamma*I)
        double beta [2] = {1.0/work->settings->gamma,0};
        work->chol->LD = cholmod_analyze (work->data->Q, &work->chol->c) ;

        cholmod_factorize_p (work->data->Q, beta, NULL, 0, work->chol->LD, &work->chol->c);
        if ((&work->chol->c)->status != CHOLMOD_OK) {
            (&work->chol->c)->supernodal = CHOLMOD_SIMPLICIAL;
            cholmod_factorize_p (work->data->Q, beta, NULL, 0, work->chol->LD, &work->chol->c);
        }

    }
    //d = ldlsolve(LD, -dphi)
    prea_vec_copy(work->dphi, work->neg_dphi, work->data->n);
    vec_mult_scalar(work->neg_dphi, -1, work->data->n);
    work->chol->d = cholmod_solve (CHOLMOD_LDLt, work->chol->LD, work->chol->neg_dphi, &work->chol->c);
    work->d = work->chol->d->x;

    //Store old active set
    prea_int_vec_copy(work->chol->active_constraints, work->chol->active_constraints_old, work->data->m);
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

void set_entering_leaving_constraints(QPALMWorkspace *work) {
    int nb_enter = 0;
    int nb_leave = 0;
    for (int i = 0; i < work->data->m; i++) {
        if (work->chol->active_constraints[i] && !work->chol->active_constraints_old[i]) {
            work->chol->enter[nb_enter] = i;
            nb_enter++;
        }
        if (!work->chol->active_constraints[i] && work->chol->active_constraints_old[i]) {
            work->chol->leave[nb_leave] = i;
            nb_leave++;
        }
    }
    work->chol->nb_enter = nb_enter;
    work->chol->nb_leave = nb_leave;
}