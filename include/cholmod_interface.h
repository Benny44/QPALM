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
void mat_vec(const cholmod_sparse *A,
             const cholmod_dense  *x,
             cholmod_dense        *y);

/* Matrix-transpose-vector multiplication
 *    y  =  A'*x 
 */
void mat_tpose_vec(const cholmod_sparse *A,
                   const cholmod_dense  *x,
                   cholmod_dense        *y);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef CHOLMOD_INTERFACE_H 