/**
 * @file lin_alg.h
 * @author Ben Hermans
 * @brief Linear algebra with vectors.
 * @details Common operations, such as vector products, infinity norm, elementwise 
 * add/product/division/max etc. are included in this file.
 */

#ifndef LIN_ALG_H
# define LIN_ALG_H


# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

# include "types.h"
#include "cholmod_interface.h"

/**
 * @name Vector functions
 * @{
 */

/**
 * Copy vector a into output.
 * @warning This function uses malloc.
 * @param a Vector
 * @param n Vector length
 * @return Copy of a
 */
c_float* vec_copy(c_float *a,
                  c_int    n);

/** 
 * Copy vector a into preallocated vector b.
 * 
 * @param a Input vector
 * @param b Output vector
 * @param n Vector length
 */
void prea_vec_copy(const c_float *a,
                   c_float       *b,
                   c_int          n);

/**
 * Copy integer vector a into preallocated vector b.
 * 
 * @param a Input vector
 * @param b Output vector
 * @param n Vector length
 */
void prea_int_vec_copy(const c_int *a,
                       c_int       *b,
                       c_int        n);

/**
 * Fill float vector with a scalar value.
 * 
 * @param a Vector
 * @param sc Value
 * @param n Vector length
 */
void vec_set_scalar(c_float *a,
                    c_float  sc,
                    c_int    n);

/**
 * Fill int vector with a scalar value.
 * 
 * @param a Vector
 * @param sc Value
 * @param n Vector length
 */
void vec_set_scalar_int(c_int *a, 
                        c_int  sc, 
                        c_int  n);

/**
 * Mulitply vector with a constant scale factor.
 * 
 * @param a Vector
 * @param sc Value
 * @param n Vector length 
 */
void vec_mult_scalar(c_float *a,
                     c_float  sc,
                     c_int    n);

/** 
 * Inner product between two vectors, @f$a^T \cdot b@f$. 
 * 
 * @param a Vector
 * @param b Vector
 * @param n Vector length
 * @return Result of the inner product
 */
c_float vec_prod(const c_float *a,
                 const c_float *b,
                 c_int          n);

/**
 * Scaled addition of one vector to another vector, @f$c_i = a_i + sc\cdot b_i@f$ 
 * 
 * @param a Input vector
 * @param b Input vector
 * @param c Output vector
 * @param sc Scaling value
 * @param n Vector length
 * */
void vec_add_scaled(const c_float *a,
                    const c_float *b,
                    c_float       *c,
                    c_float        sc,
                    c_int          n);

/**
 * Infinity norm of a vector, @f$\|a\|_\infty@f$
 * 
 * @param a Vector
 * @param n Vector length
 * @return Infinity norm of a
 */
c_float vec_norm_inf(const c_float *a,
                     c_int          n);

/**
 * Elementwise reciprocal @f$b_i = 1/a_i@f$.
 * 
 * This function is used in scaling.
 * 
 * @param a Input vector
 * @param b Output vector
 * @param n Vector length
 */
void vec_ew_recipr(const c_float *a,
                   c_float       *b,
                   c_int          n);

/** 
 * Elementwise maximum between vectors, @f$c_i = \textrm{max}(a_i, b_i)@f$.
 * 
 * @param a Input vector
 * @param b Input vector
 * @param c Output vector
 * @param n Vector length
 */
void vec_ew_max_vec(const c_float *a,
                    const c_float *b,
                    c_float       *c,
                    c_int          n);

/** 
 * Elementwise minimum between vectors, @f$c_i = \textrm{min}(a_i, b_i)@f$.
 * 
 * @param a Input vector
 * @param b Input vector
 * @param c Output vector
 * @param n Vector length
 */
void vec_ew_min_vec(const c_float *a,
                    const c_float *b,
                    c_float       *c,
                    c_int          n);

/** 
 * Elementwise mid between vectors, @f$c_i = \textrm{max}(b_{\textrm{min},i}, \textrm{min}(a_i, b_{\textrm{max},i}))@f$.
 * 
 * @param a Input vector
 * @param bmin Lower bounds
 * @param bmax Upper bounds
 * @param c Output vector
 * @param n Vector length
 */
void vec_ew_mid_vec(const c_float *a,
                    const c_float *bmin,
                    const c_float *bmax,
                    c_float       *c,
                    c_int          n);

/**
 * Elementwise product, @f$c_i = a_i\cdot b_i@f$.
 * 
 * @param a Input vector
 * @param b Input vector
 * @param c Output vector
 * @param n Vector length
 */
void vec_ew_prod(const c_float *a,
                 const c_float *b,
                 c_float       *c,
                 c_int          n);

/**
 * Elementwise division, @f$c_i = a_i/b_i@f$.
 * 
 * @param a Input vector
 * @param b Input vector
 * @param c Output vector
 * @param n Vector length
 */
void vec_ew_div(const c_float *a,
                const c_float *b,
                c_float       *c,
                c_int          n);

/** 
 * Elementwise square root, @f$b_i = \sqrt{a_i}@f$ 
 * 
 * @param a Input vector
 * @param b Output vector
 * @param n Vector length
 */
void vec_ew_sqrt(const c_float *a,
                 c_float       *b,
                 c_int          n);


/**
 * @}
 */

/* MATRIX FUNCTIONS ----------------------------------------------------------*/

/* Moved to cholmod_interface.c*/


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef LIN_ALG_H