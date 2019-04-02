#ifndef CHOLMOD_INTERFACE_H
#define CHOLMOD_INTERFACE_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "cholmod.h"
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




# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CHOLMOD_INTERFACE_H 