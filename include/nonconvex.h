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

/**
 * Calculate the minimum eigenvalue of @f$Q@f$
 * 
 * The routine can be broken down into two steps.
 * 1) Compute @f$\lambda =@f$ the maximum (in absolute value) eigenvalue of @f$Q@f$. If this is negative, then return.
 * 2) Compute the maximum (in absolute value) eigenvalue of @f$Q - \lambda*I@f$. This is the smallest eigenvalue
 * of @f$Q@f$ offset by the initial @f$\lambda@f$, so return after adjusting.
 * 
 * @param work Workspace
 */
c_float minimum_eigenvalue_Q(QPALMWorkspace *work);

/**
 * Calculate the maximum (in absolute value) eigenvalue of @f$ Q + \gamma*I @f$
 * 
 * This routine implements a simple version of power iterations. 
 * 
 * @param work Workspace
 * @param gamma Offset 
 */
c_float power_iterations_Q(QPALMWorkspace *work, 
                         c_float         gamma);

/**
 * Calculate the Gershgorin upper bound for the eigenvalues of a symmetric matrix.
 * 
 * This routine uses the Gershgorin circle theorem to compute an upper bound on the eigenvalues of 
 * a matrix.
 * 
 * @assumption M is symmetric
 * @param M Matrix
 * @param center Vector of size M->ncol to hold the values of the centers of the discs
 * @param radius Vector of size M->ncol to hold the values of the radii of the discs
 * @return Upper bound on the eigenvalues of M
 */
c_float gershgorin_max(cholmod_sparse* M, c_float *center, c_float *radius);

/**
 * Set the proximal parameters for nonconvex QPs.
 * 
 * QPALM can deal with nonconvex QPs, by setting the initial and maximal proximal penalty small enough 
 * (smaller than @f$ \frac{1}{|\lambda_\textrm{min}|} @f$). This ensures positive definiteness of 
 * @f$ Q + \frac{1}{\gamma}I @f$ during the iterations.
 * 
 * @param work Workspace
 */
void set_settings_nonconvex(QPALMWorkspace *work);


# ifdef __cplusplus
}
# endif /* ifdef __cplusplus */

#endif /*NONCONVEX_H*/