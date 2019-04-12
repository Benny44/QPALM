/**
 * @file newton.c
 * @author Ben Hermans
 * @brief Functions to calculate the semismooth Newton direction.
 * @details The functions in this file concern the calculation of the semismooth Newton direction. 
 * Factorizing, updating the factorization and solving the linear system are performed by functions in 
 * cholmod_interface.c. 
 */

#include "newton.h"
#include "lin_alg.h"
#include "cholmod.h"
#include <stdio.h>

void newton_set_direction(QPALMWorkspace *work) {

    set_active_constraints(work);
    if (work->chol->reset_newton && work->chol->nb_active_constraints) {
        work->chol->reset_newton = FALSE;
        ldlcholQAtsigmaA(work);   
    } else if (work->chol->nb_active_constraints) {
        set_entering_leaving_constraints(work);
        if(work->chol->nb_enter) {
            ldlupdate_entering_constraints(work);
        }
        if(work->chol->nb_leave) {
            ldldowndate_leaving_constraints(work); 
        }
    } else {
        ldlcholQ(work);
    }
    ldlsolveLD_neg_dphi(work);

    //Store old active set
    prea_int_vec_copy(work->chol->active_constraints, work->chol->active_constraints_old, work->data->m);
}

void set_active_constraints(QPALMWorkspace *work) {
    work->chol->nb_active_constraints = 0;
    for (size_t i = 0; i < work->data->m; i++) {
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
    for (size_t i = 0; i < work->data->m; i++) {
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