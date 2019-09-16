/**
 * @file lin_alg.c
 * @author Ben Hermans
 * @brief Linear algebra with vectors.
 * @details Common operations, such as vector products, infinity norm, elementwise 
 * add/product/division/max etc. are included in this file.
 */

#include "lin_alg.h"

c_float* vec_copy(const c_float *a, size_t n) {
  c_float *b;
  size_t    i;

  b = c_malloc(n * sizeof(c_float));

  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }

  return b;
}

void prea_vec_copy(const c_float *a, c_float *b, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

void prea_int_vec_copy(const c_int *a, c_int *b, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    b[i] = a[i];
  }
}

void vec_set_scalar(c_float *a, c_float sc, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    a[i] = sc;
  }
}

void vec_set_scalar_int(c_int *a, c_int sc, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    a[i] = sc;
  }
}

void vec_mult_scalar(c_float *a, c_float sc, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    a[i] *= sc;
  }
}

c_float vec_prod(const c_float *a, const c_float *b, size_t n) {
  c_float prod = 0.0;
  size_t i = 0; 

  if(n >= 4) {
      for (; i <= n-4; i+=4) {
        prod += (a[i]*b[i] + a[i+1]*b[i+1] + a[i+2]*b[i+2] + a[i+3]*b[i+3]);
      }
  }
  for (; i < n; i++) {
    prod += a[i] * b[i];
  }

  return prod;
}


void vec_ew_prod(const c_float *a, const c_float *b, c_float *c, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    c[i] = a[i] * b[i];
  }
}


void vec_ew_div(const c_float *a, const c_float *b, c_float *c, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    c[i] = a[i] / b[i];
  }
}


void vec_add_scaled(const c_float *a, const c_float *b, c_float *c, c_float sc, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    c[i] =  a[i] + sc * b[i];
  }
}

// c_float vec_norm_inf(const c_float *a, size_t n) {
//   size_t   i;
//   c_float abs_a_i;
//   c_float max = 0.0;

//   for (i = 0; i < n; i++) {
//     abs_a_i = c_absval(a[i]);

//     if (abs_a_i > max) max = abs_a_i;
//   }
//   return max;
// }

c_float vec_norm_inf(const c_float *a, size_t n) {
    register size_t j = 0;
    register c_float s0 = 0.;
    register c_float s1 = 0.;
    register c_float s2 = 0.;
    register c_float s3 = 0.;
    register c_float max0 = 0.;
    register c_float max1 = 0.;
    register c_float max2 = 0.;
    register c_float max3 = 0.;
    const size_t block_size = 4;
    const size_t block_len = n >> 2;
    const size_t remaining = n % block_size; /*Initializing four blocks for 
                                                * efficient implementation of 
                                                * inner product using registers */

    while (j < block_len * block_size) {
      s0 = c_absval(a[j]); max0 = s0 > max0 ? s0 : max0;
      s1 = c_absval(a[j+1]); max1 = s1 > max1 ? s1 : max1;
      s2 = c_absval(a[j+2]); max2 = s2 > max2 ? s2 : max2;
      s3 = c_absval(a[j+3]); max3 = s3 > max3 ? s3 : max3;
      j+=4;
    }    

    max0 = max0 > max1 ? max0 : max1;
    max0 = max0 > max2 ? max0 : max2;
    max0 = max0 > max3 ? max0 : max3;
    j = block_size * block_len;
    switch (remaining) {
        case 3: max0 = max0 > c_absval(a[j+2]) ? max0 : c_absval(a[j+2]);
        case 2: max0 = max0 > c_absval(a[j+1]) ? max0 : c_absval(a[j+1]);
        case 1: max0 = max0 > c_absval(a[j+0]) ? max0 : c_absval(a[j]); /*Taking contribution from the last terms
                                    * that were not included in the block*/
        case 0:;
    }
    return max0;

}

void vec_ew_recipr(const c_float *a, c_float *b, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    b[i] = (c_float)1.0 / a[i];
  }
}

void vec_ew_max_vec(const c_float *a, const c_float *b, c_float *c, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    c[i] = c_max(a[i], b[i]);
  }
}

void vec_ew_min_vec(const c_float *a, const c_float *b, c_float *c, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    c[i] = c_min(a[i], b[i]);
  }
}

void vec_ew_mid_vec(const c_float *a, const c_float *bmin, const c_float *bmax, c_float *c, size_t n) {
  size_t i;
  
  for (i = 0; i < n; i++) {
    c[i] = c_max(bmin[i], c_min(a[i], bmax[i]));
  }
}

void vec_ew_sqrt(const c_float *a, c_float *b, size_t n) {
  size_t i;

  for (i = 0; i < n; i++) {
    b[i] = c_sqrt(a[i]);
  }
}

/**
 * @}
 */


/* MATRIX FUNCTIONS ----------------------------------------------------------*/

/* Moved to cholmod_interface.c*/

