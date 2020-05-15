/**
 * @file newton.c
 * @author Ben Hermans
 * @brief Functions to calculate the semismooth Newton direction.
 * @details The functions in this file concern the calculation of the semismooth Newton direction. 
 * Factorizing, updating the factorization and solving the linear system are performed by functions in 
 * cholmod_interface.c. 
 */
# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "newton.h"
#include "lin_alg.h"
#include <stdio.h>

#ifdef USE_CHOLMOD
#include "cholmod.h"
void newton_set_direction(QPALMWorkspace *work, solver_common *c) {

    set_active_constraints(work);
    set_entering_leaving_constraints(work);
    if ((work->solver->reset_newton && work->solver->nb_active_constraints) || 
        (work->solver->nb_enter + work->solver->nb_leave) > MAX_RANK_UPDATE) {
        work->solver->reset_newton = FALSE;
        ldlcholQAtsigmaA(work, c);   
    } else if (work->solver->nb_active_constraints) {
        if(work->solver->nb_enter) {
            ldlupdate_entering_constraints(work, c);
        }
        if(work->solver->nb_leave) {
            ldldowndate_leaving_constraints(work, c); 
        }
    } else {
        ldlchol(work->data->Q, work, c);
    }

    ldlsolveLD_neg_dphi(work, c);

    //Store old active set
    prea_int_vec_copy(work->solver->active_constraints, work->solver->active_constraints_old, work->data->m);
}

#endif /* USE_CHOLMOD */

void set_active_constraints(QPALMWorkspace *work) {
    work->solver->nb_active_constraints = 0;
    for (size_t i = 0; i < work->data->m; i++) {
        if ((work->Axys[i] <= work->data->bmin[i]) || ((work->Axys[i] >= work->data->bmax[i]))){
            work->solver->active_constraints[i] = TRUE;
            work->solver->nb_active_constraints++;
        } else {
            work->solver->active_constraints[i] = FALSE;
        }         
    }
}

void set_entering_leaving_constraints(QPALMWorkspace *work) {
    int nb_enter = 0;
    int nb_leave = 0;
    for (size_t i = 0; i < work->data->m; i++) {
        if (work->solver->active_constraints[i] && !work->solver->active_constraints_old[i]) {
            work->solver->enter[nb_enter] = (c_int)i;
            nb_enter++;
        }
        if (!work->solver->active_constraints[i] && work->solver->active_constraints_old[i]) {
            work->solver->leave[nb_leave] = (c_int)i;
            nb_leave++;
        }
    }
    work->solver->nb_enter = nb_enter;
    work->solver->nb_leave = nb_leave;
}



# ifdef __cplusplus
}
# endif // ifdef __cplusplus