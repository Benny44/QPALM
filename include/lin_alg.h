#ifndef LIN_ALG_H
# define LIN_ALG_H


# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

# include "types.h"


/* VECTOR FUNCTIONS ----------------------------------------------------------*/

/* copy vector a into output (Uses MALLOC)*/
c_float* vec_copy(c_float *a,
                  c_int    n);

/* copy vector a into preallocated vector b */
void prea_vec_copy(const c_float *a,
                   c_float       *b,
                   c_int          n);

/* copy integer vector a into preallocated vector b */
void prea_int_vec_copy(const c_int *a,
                       c_int       *b,
                       c_int        n);

/* set float vector to scalar */
void vec_set_scalar(c_float *a,
                    c_float  sc,
                    c_int    n);

/* multiply scalar to vector */
void vec_mult_scalar(c_float *a,
                     c_float  sc,
                     c_int    n);

/* Inner product a'b */
c_float vec_prod(const c_float *a,
                 const c_float *b,
                 c_int          n);

/* c = a + sc*b */
void vec_add_scaled(const c_float *a,
                    const c_float *b,
                    c_float       *c,
                    c_float        sc,
                    c_int          n);

/* ||v||_inf */
c_float vec_norm_inf(const c_float *v,
                     c_int          l);

/* Vector elementwise reciprocal b = 1./a (needed for scaling)*/
void vec_ew_recipr(const c_float *a,
                   c_float       *b,
                   c_int          n);

/* Elementwise maximum between vectors, c = max(a, b) */
void vec_ew_max_vec(const c_float *a,
                    const c_float *b,
                    c_float       *c,
                    c_int          n);

/* Elementwise minimum between vectors, c = min(a, b) */
void vec_ew_min_vec(const c_float *a,
                    const c_float *b,
                    c_float       *c,
                    c_int          n);

/* Elementwise mid between vectors, c = max(bmin, min(a, bmax)) */
void vec_ew_mid_vec(const c_float *a,
                    const c_float *bmin,
                    const c_float *bmax,
                    c_float       *c,
                    c_int          n);

/* Elementwise product a.*b stored in c*/
void vec_ew_prod(const c_float *a,
                 const c_float *b,
                 c_float       *c,
                 c_int          n);

/* Elementwise division a./b stored in c*/
void vec_ew_div(const c_float *a,
                const c_float *b,
                c_float       *c,
                c_int          n);

/* Elementwise sqrt of the vector elements b = sqrt(a) */
void vec_ew_sqrt(const c_float *a,
                 c_float       *b,
                 c_int          n);



/* MATRIX FUNCTIONS ----------------------------------------------------------*/


/* Premultiply matrix A by diagonal matrix with diagonal d,
   i.e. scale the rows of A by d
 */
void mat_premult_diag(csc           *A,
                      const c_float *d);

/* Premultiply matrix A by diagonal matrix with diagonal d,
   i.e. scale the columns of A by d
 */
void mat_postmult_diag(csc           *A,
                       const c_float *d);

/* Matrix-vector multiplication
 *    y  =  A*x 
 */
void mat_vec(const csc     *A,
             const c_float *x,
             c_float       *y);


/* Matrix-transpose-vector multiplication
 *    y  =  A'*x 
 */
void mat_tpose_vec(const csc     *A,
                   const c_float *x,
                   c_float       *y);

/* Matrix-transpose-vector multiplication while skipping the diagonal
 *    y +=  A'*x  (if plus_eq == 1)
 */
void mat_tpose_vec_skip_diag_add(const csc     *A,
                                 const c_float *x,
                                 c_float       *y);

/* Matrix-vector multiplication (where only Q's upper triagonal elements are stored)
 *    y  =  Q*x  
 */
void mat_vec_triu(const csc     *Q,
                  const c_float *x,
                  c_float       *y);

/**
 * Infinity norm of each matrix column
 * @param M	Input matrix
 * @param E     Vector of infinity norms
 *
 */
void mat_inf_norm_cols(const csc *M,
                       c_float   *E);

/**
 * Infinity norm of each matrix row
 * @param M	Input matrix
 * @param E     Vector of infinity norms
 *
 */
void mat_inf_norm_rows(const csc *M,
                       c_float   *E);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef LIN_ALG_H