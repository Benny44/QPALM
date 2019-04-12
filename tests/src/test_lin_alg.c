#include "lin_alg.h"
#include "global_opts.h"
#include <CUnit/CUnit.h>

#define DIM 3
#define TOL 1e-8

c_float a[DIM] = {0.1, 2.5, -3.9};
c_float b[DIM] = {0.0, 10, 4};
c_float c[DIM] = {0.0, 0.0, 0.0};

void reset_abc(void) {
    a[0] = 0.1; a[1] = 2.5; a[2] = -3.9;
    b[0] = 0.0; b[1] = 10; b[2] = 4;
    c[0] = 0.0; c[1] = 0.0; c[2] = 0.0;
}

void test_vec_set_scalar(void) {
    c_float scalar = 5.5;
    vec_set_scalar(a, scalar, DIM);
    for (size_t i = 0; i < DIM; i++) {
        CU_ASSERT_DOUBLE_EQUAL(a[i], scalar, TOL);
    }
    reset_abc();
}

void test_vec_set_scalar_int(void){
    c_int scalar = 3;
    vec_set_scalar(a, scalar, DIM);
    for (size_t i = 0; i < DIM; i++) {
        CU_ASSERT_EQUAL(a[i], scalar);
    }
    reset_abc();
}

void test_vec_mult_scalar(void);

void test_vec_prod(void);

void test_vec_add_scaled(void);

void test_vec_norm_inf(void);

void test_vec_ew_recipr(void);

void test_vec_ew_max_vec(void);

void test_vec_ew_min_vec(void);

void test_vec_ew_mid_vec(void);

void test_vec_ew_prod(void);

void test_vec_ew_div(void);

void test_vec_ew_sqrt(void);