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

#ifdef USE_LADEL
#include "ladel.h"
#include "ladel_copy.h"
#include "ladel_global.h"
#include "ladel_types.h"
#include "ladel_row_mod.h"
#include "ladel_debug_print.h"
#elif defined USE_CHOLMOD
#include "cholmod.h"
#endif

#ifdef USE_LADEL
static void qpalm_form_kkt(QPALMWorkspace *work)
{
    solver_sparse *Q = work->data->Q, *A = work->data->A, *kkt = work->solver->kkt, *kkt_full = work->solver->kkt_full, *At = work->solver->At;
    ladel_int col, index, index_kkt, n = work->data->n, m = work->data->m, Qnz = Q->nzmax;
    c_float *sigma_inv = work->sigma_inv, *first_elem_A = work->solver->first_elem_A;
    c_int *first_row_A = work->solver->first_row_A;
    /* copy Q */
    for (col = 0; col < n; col++)
    {
        kkt->p[col] = kkt_full->p[col] = Q->p[col];
        kkt->nz[col] = Q->p[col+1] - Q->p[col];
    }
    kkt->p[col] = kkt_full->p[col] = Q->p[col];
    for (index = 0; index < Qnz; index++)
    {
        kkt->i[index] = kkt_full->i[index] = Q->i[index];
        kkt->x[index] = kkt_full->x[index] = Q->x[index];
    }

    /* copy [At; -\Sigma^{-1}] */
    index_kkt = Qnz;
    for (; col < m+n; col++)
    {
        kkt_full->i[index_kkt] = first_row_A[col-n] = At->i[At->p[col-n]];
        kkt_full->x[index_kkt] = first_elem_A[col-n] = At->x[At->p[col-n]];

        if (work->solver->active_constraints[col-n])
        {
            kkt->nz[col] = At->p[col-n+1] - At->p[col-n] + 1;
            kkt->i[index_kkt] = At->i[At->p[col-n]];
            kkt->x[index_kkt] = At->x[At->p[col-n]];
        } 
        else 
        {
            kkt->nz[col] = 1;
            kkt->i[index_kkt] = col;
            kkt->x[index_kkt] = 1;
        }

        if (At->p[col-n+1]-At->p[col-n] != 0) index_kkt++;

        for (index = At->p[col-n]+1; index < At->p[col-n+1]; index++)
        {
            kkt->i[index_kkt] = kkt_full->i[index_kkt] = At->i[index];
            kkt->x[index_kkt] = kkt_full->x[index_kkt] = At->x[index];
            index_kkt++;
        }

        kkt->i[index_kkt] = kkt_full->i[index_kkt] = col;
        kkt->x[index_kkt] = kkt_full->x[index_kkt] = -sigma_inv[col-n];
        if (At->p[col-n+1]-At->p[col-n] == 0) kkt->x[index_kkt] = 1;
        index_kkt++;

        kkt->p[col+1] = kkt_full->p[col+1] = Qnz + At->p[col+1-n] + 1 + col - n;
    }
}


static void qpalm_reform_kkt(QPALMWorkspace *work)
{
    solver_sparse *kkt = work->solver->kkt, *At = work->solver->At;
    ladel_int col, index, n = work->data->n, m = work->data->m;
    ladel_int *first_row_A = work->solver->first_row_A;
    ladel_double *sigma_inv = work->sigma_inv, *first_elem_A = work->solver->first_elem_A;

    for (col = n; col < n+m; col++)
    {
        if (work->solver->active_constraints[col-n])
        {
            kkt->nz[col] = At->p[col-n+1] - At->p[col-n] + 1;
            kkt->i[kkt->p[col]] = first_row_A[col-n];
            kkt->x[kkt->p[col]] = first_elem_A[col-n];
            kkt->x[kkt->p[col+1]-1] = -sigma_inv[col-n]; /* row index should be correct already */
            kkt->i[kkt->p[col+1]-1] = col;
        } else
        {
            kkt->nz[col] = 1;
            kkt->i[kkt->p[col]] = col;
            kkt->x[kkt->p[col]] = 1; 
        }
    }
}

static void kkt_update_entering_constraints(QPALMWorkspace *work, solver_common *c)
{
    solver_sparse *kkt = work->solver->kkt, *At = work->solver->At;
    ladel_int col, index, n = work->data->n, m = work->data->m;
    ladel_int *first_row_A = work->solver->first_row_A;
    ladel_double *sigma_inv = work->sigma_inv, *first_elem_A = work->solver->first_elem_A;

    for (index = 0; index < work->solver->nb_enter; index++)
    {
        col = work->solver->enter[index] + n;
        kkt->nz[col] = At->p[col-n+1] - At->p[col-n] + 1;
        kkt->i[kkt->p[col]] = first_row_A[col-n];
        kkt->x[kkt->p[col]] = first_elem_A[col-n];
        kkt->x[kkt->p[col+1]-1] = -sigma_inv[col-n]; /* row index should be correct already */
        ladel_row_add(work->solver->LD, work->solver->sym, col, kkt, col, -sigma_inv[col-n], c);
    }
}

static void kkt_update_leaving_constraints(QPALMWorkspace *work, solver_common *c)
{
    ladel_int col, index, n = work->data->n, m = work->data->m;

    for (index = 0; index < work->solver->nb_leave; index++)
    {
        col = work->solver->leave[index] + n;
        ladel_row_del(work->solver->LD, work->solver->sym, col, c);
        /* no need to update the kkt system here */
    }
}

static void kkt_solve(QPALMWorkspace *work, solver_common *c)
{
    c_int n = work->data->n, m = work->data->m;
    prea_vec_copy(work->dphi, work->solver->rhs_kkt, n);
    vec_self_mult_scalar(work->solver->rhs_kkt, -1, n);
    vec_set_scalar(work->solver->rhs_kkt + n, 0, m);

    ladel_dense_solve(work->solver->LD, work->solver->rhs_kkt, work->solver->sol_kkt, c);
    prea_vec_copy(work->solver->sol_kkt, work->d, n);
}

#endif

void newton_set_direction(QPALMWorkspace *work, solver_common *c) {

    set_active_constraints(work);
    set_entering_leaving_constraints(work);
    #ifdef USE_LADEL
    ladel_diag d;
    d.diag_elem = 1.0/work->gamma;
    if (work->settings->proximal) d.diag_size = work->data->n;
    else d.diag_size = 0;

    if (work->solver->first_factorization)
    {
        qpalm_form_kkt(work);
        ladel_factorize_advanced_with_diag(work->solver->kkt, d, work->solver->sym, work->settings->ordering, &work->solver->LD, work->solver->kkt_full, c);
        work->solver->first_factorization = FALSE;
    } 
    else if (work->solver->reset_newton || 
            (work->solver->nb_enter + work->solver->nb_leave) > 0.1*(work->data->n+work->data->m)) 
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

    #elif defined USE_CHOLMOD
    if ((work->solver->reset_newton && work->solver->nb_active_constraints) || 
        (work->solver->nb_enter + work->solver->nb_leave) > MAX_RANK_UPDATE) {
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
    #endif /* USE_CHOLMOD */

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