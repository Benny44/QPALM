#ifndef CHOLMOD_INTERFACE_H
#define CHOLMOD_INTERFACE_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "cholmod.h"
#include "global_opts.h"
#include "types.h"


/* Matrix-vector multiplication
 *    y  =  A*x 
 */
void mat_vec(cholmod_sparse *A,
             cholmod_dense  *x,
             cholmod_dense        *y,
             cholmod_common       *c);

/* Matrix-transpose-vector multiplication
 *    y  =  A'*x 
 */
void mat_tpose_vec(cholmod_sparse *A,
                   cholmod_dense  *x,
                   cholmod_dense        *y,
                   cholmod_common       *c);

/**
 * Infinity norm of each matrix column
 * @param M	Input matrix
 * @param E     Vector of infinity norms
 *
 */
void mat_inf_norm_cols(cholmod_sparse *M,
                       c_float   *E);

/**
 * Infinity norm of each matrix row
 * @param M	Input matrix
 * @param E     Vector of infinity norms
 *
 */
void mat_inf_norm_rows(cholmod_sparse *M,
                       c_float   *E);

/**
 * Calculate LDL factorization of M
 * If work->settings->proximal = true, use M+(1/gamma)*I instead
 * @param M matrix to be factorized
 * @param work qpalm workspace
 */
void ldlchol(cholmod_sparse *M, 
             QPALMWorkspace *work);

/**
 * Calculate LDL factorization of Q
 * If work->settings->proximal = true, use Q+(1/gamma)*I instead
 * @param work qpalm workspace
 */
void ldlcholQ(QPALMWorkspace *work);

/**
 * Calculate LDL factorization of Q+A(active_cnstrs,:)'*diag(sigma)*A(active_cnstrs,:)
 * If work->settings->proximal = true, use Q+(1/gamma)*I+A(active_cnstrs,:)'*diag(sigma)*A(active_cnstrs,:) instead
 * @param work qpalm workspace
 */
void ldlcholQAtsigmaA(QPALMWorkspace *work);

/**
 * Update the LDL factorization given entering constraints
 * @param work qpalm workspace
 */
void ldlupdate_entering_constraints(QPALMWorkspace *work);

/**
 * Downdate the LDL factorization given leaving constraints
 * @param work qpalm workspace
 */
void ldldowndate_leaving_constraints(QPALMWorkspace *work);

/**
 * Solve the linear system LDL'd = -dphi
 * @param work qpalm workspace
 */
void ldlsolveLD_neg_dphi(QPALMWorkspace *work);

/**
 * Cholmod settings to indicate that we use an LDLt factorization
 * @param c Cholmod workspace
 */
void cholmod_set_settings(cholmod_common *c);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CHOLMOD_INTERFACE_H 