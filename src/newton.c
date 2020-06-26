/**
 * @file newton.c
 * @author Ben Hermans
 * @brief Functions to calculate the semismooth Newton direction.
 * @details The functions in this file concern the calculation of the semismooth Newton direction. 
 * Factorizing, updating the factorization and solving the linear system are performed by functions in 
 * solver_interface.c. 
 */
# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "newton.h"
#include "lin_alg.h"
#include <stdio.h>

void newton_set_direction(QPALMWorkspace *work, solver_common *c) {

    set_active_constraints(work);
    set_entering_leaving_constraints(work);
    if (work->solver->factorization_method == FACTORIZE_KKT)
    {
        #ifdef USE_LADEL 
        // NB FACTORIZE_KKT only defined with LADEL for now
        // TODO (optionally) extract the ladel functions and bury them, use cholmod also for kkt
        ladel_diag d;
        d.diag_elem = 1.0/work->gamma;
        if (work->settings->proximal) d.diag_size = work->data->n;
        else d.diag_size = 0;

        if (work->solver->first_factorization)
        {
            qpalm_form_kkt(work);
            work->solver->LD = ladel_factor_free(work->solver->LD);
            ladel_factorize_advanced_with_diag(work->solver->kkt, d, work->solver->sym, work->settings->ordering, &work->solver->LD, work->solver->kkt_full, c);
            work->solver->first_factorization = FALSE;
        } 
        else if (work->solver->reset_newton || 
                (work->solver->nb_enter + work->solver->nb_leave) > c_min(0.1*(work->data->n+work->data->m), MAX_RANK_UPDATE)) 
        {
            qpalm_reform_kkt(work);
            ladel_factorize_with_prior_basis_with_diag(work->solver->kkt, d, work->solver->sym, work->solver->LD, c);
        }
        else 
        {
            if(work->solver->nb_enter) 
                kkt_update_entering_constraints(work, c);

            if (work->solver->nb_leave)
                kkt_update_leaving_constraints(work, c);
        }

        kkt_solve(work, c);
        #endif
    } else if (work->solver->factorization_method == FACTORIZE_SCHUR)
    {
        if ((work->solver->reset_newton && work->solver->nb_active_constraints) || 
            (work->solver->nb_enter + work->solver->nb_leave) > c_min(0.1*(work->data->n+work->data->m), MAX_RANK_UPDATE)) {
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
    } 
    
    //Store old active set
    prea_int_vec_copy(work->solver->active_constraints, work->solver->active_constraints_old, work->data->m);

    work->solver->reset_newton = FALSE;
}

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