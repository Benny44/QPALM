/**
 * @file cholmod_interface.h
 * @author Ben Hermans
 * @brief Interface and wrapper to cholmod functions
 * @details This file includes all calls to cholmod functions apart from scaling in scaling.c and memory
 * allocation/deallocation in the main functions in qpalm.c. It includes all matrix operations, such as
 * matrix vector products, row- and columnwise norms, cholesky factorizations, factorization updates and
 * solving the linear system. Finally, all the settings relevant to cholmod (and suitesparse) are included
 * in this file as well.
 */

#ifndef CHOLMOD_INTERFACE_H
#define CHOLMOD_INTERFACE_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "cholmod.h"
#include "global_opts.h"
#include "constants.h"
#include "types.h"


/** 
 * Matrix-vector multiplication.
 *    
 * y  =  A*x
 * 
 * @param A Sparse matrix
 * @param x Dense input vector
 * @param y Dense output vector
 * @param c Cholmod environment 
 */
void mat_vec(cholmod_sparse *A,
             cholmod_dense  *x,
             cholmod_dense        *y,
             cholmod_common       *c);

/** 
 * Matrix-transpose-vector multiplication.
 *    
 * y  =  A'*x 
 * 
 * @param A Sparse matrix
 * @param x Dense input vector
 * @param y Dense output vector
 * @param c Cholmod environment 
 */
void mat_tpose_vec(cholmod_sparse *A,
                   cholmod_dense  *x,
                   cholmod_dense        *y,
                   cholmod_common       *c);

/**
 * Infinity norm of each matrix column.
 * 
 * @param M	Sparse input matrix
 * @param E Vector of infinity norms
 *
 */
void mat_inf_norm_cols(cholmod_sparse *M,
                       c_float   *E);

/**
 * Infinity norm of each matrix row.
 * 
 * @param M	Sparse input matrix
 * @param E Vector of infinity norms
 *
 */
void mat_inf_norm_rows(cholmod_sparse *M,
                       c_float   *E);

/**
 * Calculate LDL factorization of a matrix M.
 * 
 * If work->settings->proximal = true, use M+(1/gamma)*I instead.
 * 
 * @param M Matrix to be factorized
 * @param work Workspace
 */
void ldlchol(cholmod_sparse *M, 
             QPALMWorkspace *work);

/**
 * Calculate LDL factorization of Q.
 * 
 * Q is the quadratic part of the objective function.
 * 
 * If work->settings->proximal = true, use Q+(1/gamma)*I instead.
 * 
 * @param work Workspace
 */
void ldlcholQ(QPALMWorkspace *work);

/**
 * Calculate LDL factorization of Q+A_a'*Sigma_{a,a}*A_a, with Sigma=diag(sigma) and a the set of active constraints.
 * 
 * If work->settings->proximal = true, use Q+(1/gamma)*I+A_a'*Sigma_{a,a}*A_a instead.
 * 
 * @param work Workspace
 */
void ldlcholQAtsigmaA(QPALMWorkspace *work);

/**
 * Update the LDL factorization given a set of entering constraints.
 * 
 * The index set of entering constraints is assumed to be set in work->chol->enter.
 * 
 * @param work Workspace
 */
void ldlupdate_entering_constraints(QPALMWorkspace *work);

/**
 * Downdate the LDL factorization given a set of leaving constraints.
 * 
 * The index set of leaving constraints is assumed to be set in work->chol->leave.
 * 
 * @param work Workspace
 */
void ldldowndate_leaving_constraints(QPALMWorkspace *work);

/**
 * Solve the linear system LDL'*d = -dphi.
 * 
 * @param work Workspace
 */
void ldlsolveLD_neg_dphi(QPALMWorkspace *work);

/**
 * Cholmod settings to indicate that we use an LDL' factorization.
 * 
 * In addition, this function sets the suitesparse memory allocation function to our custom functions.
 * 
 * @param c Cholmod environment
 */
void cholmod_set_settings(cholmod_common *c);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CHOLMOD_INTERFACE_H 