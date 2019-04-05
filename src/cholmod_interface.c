#include "cholmod_interface.h"
#include <stdio.h>

void mat_vec(cholmod_sparse *A, cholmod_dense *x, cholmod_dense *y, cholmod_common *c) {
    double one [2] = {1,0};
    double zero [2] = {0,0};
    cholmod_sdmult(A, 0, one, zero, x, y, c);
}

void mat_tpose_vec(cholmod_sparse *A, cholmod_dense *x, cholmod_dense *y, cholmod_common *c) {
    double one [2] = {1,0};
    double zero [2] = {0,0};
    cholmod_sdmult(A, 1, one, zero, x, y, c);
}

void mat_inf_norm_cols(cholmod_sparse *M, c_float *E) {
  int j, i;
  c_float *Mx = M->x;
  int *Mp = M->p;

  // Initialize zero max elements
  for (j = 0; j < M->ncol; j++) {
    E[j] = 0.;
  }

  // Compute maximum across columns
  for (j = 0; j < M->ncol; j++) {
    for (i = Mp[j]; i < Mp[j + 1]; i++) {
      E[j] = c_max(c_absval(Mx[i]), E[j]);
    }
  }
}

void mat_inf_norm_rows(cholmod_sparse *M, c_float *E) {
  int i, j, k;
  int *Mp = M->p;
  int *Mi = M->i;
  c_float *Mx = M->x;
 
  // Initialize zero max elements
  for (j = 0; j < M->nrow; j++) {
    E[j] = 0.;
  }

  // Compute maximum across rows
  for (j = 0; j < M->ncol; j++) {
    for (k = Mp[j]; k < Mp[j + 1]; k++) {
      i    = Mi[k];
      E[i] = c_max(c_absval(Mx[k]), E[i]);
    }
  }
}

void ldlchol(cholmod_sparse *M, QPALMWorkspace *work) {
  work->chol->LD = cholmod_analyze (M, &work->chol->c) ;
  if (work->settings->proximal) {
    double beta [2] = {1.0/work->settings->gamma,0};
    cholmod_factorize_p (M, beta, NULL, 0, work->chol->LD, &work->chol->c);
  } else {
    cholmod_factorize (M, work->chol->LD, &work->chol->c);
  }
  if ((&work->chol->c)->status != CHOLMOD_OK) {
      (&work->chol->c)->supernodal = CHOLMOD_SIMPLICIAL;
      if (work->settings->proximal) {
        double beta [2] = {1.0/work->settings->gamma,0};
        cholmod_factorize_p (M, beta, NULL, 0, work->chol->LD, &work->chol->c);
      } else {
        cholmod_factorize (M, work->chol->LD, &work->chol->c);
      }
  }
}

void ldlcholQ(QPALMWorkspace *work) {
  ldlchol(work->data->Q, work);
}

void ldlcholQAtsigmaA(QPALMWorkspace *work) {
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
  QAtsigmaA->stype = work->data->Q->stype;
  
  ldlchol(QAtsigmaA, work);

  cholmod_free_sparse(&AtsigmaA, &work->chol->c);
  cholmod_free_sparse(&QAtsigmaA, &work->chol->c);
}

void ldlupdate_entering_constraints(QPALMWorkspace *work) {
  cholmod_sparse *Ae;
  Ae = cholmod_submatrix(work->chol->At_sqrt_sigma, NULL, -1, 
                      work->chol->enter, work->chol->nb_enter, TRUE, TRUE, &work->chol->c);
  //LD = ldlupdate(LD,Ae,'+');
  cholmod_updown(TRUE, Ae, work->chol->LD, &work->chol->c);
  cholmod_free_sparse(&Ae, &work->chol->c);
}

void ldldowndate_leaving_constraints(QPALMWorkspace *work) {
  cholmod_sparse *Al;
  Al = cholmod_submatrix(work->chol->At_sqrt_sigma, NULL, -1, 
                      work->chol->leave, work->chol->nb_leave, TRUE, TRUE, &work->chol->c);
  //LD = ldlupdate(LD,Ae,'+');
  cholmod_updown(FALSE, Al, work->chol->LD, &work->chol->c);
  cholmod_free_sparse(&Al, &work->chol->c);
}



void cholmod_set_settings(cholmod_common *c) {
  c->final_asis = FALSE ;
  c->final_super = FALSE ;
  c->final_ll = FALSE ;
  c->final_pack = TRUE ;
  c->final_monotonic = TRUE ;
  c->final_resymbol = TRUE ;
  c->quick_return_if_not_posdef = TRUE;
  c->nmethods = 1 ;
	c->method [0].ordering = CHOLMOD_NATURAL ;
	c->postorder = FALSE ;
}