#ifndef LINESEARCH_H
#define LINESEARCH_H

#include "types.h"
#include "constants.h"


/**
 * Execute exact linesearch (using qsort)
 * 
 * @param work Workspace
 * @return tau Stepsize
 */
c_float exact_linesearch(QPALMWorkspace *work);

/**
 * Execute exact linesearch (using Newton's method)
 * 
 * @param work Workspace
 * @return tau Stepsize
 */
c_float exact_linesearch_newton(QPALMWorkspace *work);

/**
 * Execute backtracking linesearch with Armijo condition
 * 
 * @param work Workspace
 * @return tau Stepsize
 */
c_float armijo_linesearch(QPALMWorkspace *work);

/**
 * Helper function to copy a vector a in an array b (with indices)
 */
void vec_array_copy(c_float       *a, 
                    array_element *b, 
                    c_int          n);

/**
 * Select subsequence based on index list (fill the rest with QPALM_INFTY)
 */
void select_subsequence(const array_element *a, 
                        array_element       *b,
                        const c_int         *L,
                        c_int                n);

/* Inner product over index set L,  a(L)'b(L) */
c_float vec_prod_ind(const c_float *a,
                     const c_float *b,
                     const c_int   *L,
                     c_int          n);

/**
 * Helper function for qsort
 * Carries out comparison between two array_elements
 */
int compare (const void * a, const void * b);

#endif