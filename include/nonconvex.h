/**
 * @file nonconvex.h
 * @author Ben Hermans
 * @brief Routines to deal with nonconvex QPs.
 * @details The functions in this file serve to set up QPALM for a nonconvex QP. The main routine in 
 * this file computes the minimum eigenvalue of a square matrix, based on power iterations. Furthermore, 
 * some setting updates are performed. 
 */

#ifndef NONCONVEX_H
#define NONCONVEX_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"

c_float minimum_eigenvalue_Q(QPALMWorkspace *work);

c_float power_iterations_Q(QPALMWorkspace *work, 
                         c_float         gamma);

void set_settings_nonconvex(QPALMWorkspace *work);


# ifdef __cplusplus
}
# endif /* ifdef __cplusplus */

#endif /*NONCONVEX_H*/