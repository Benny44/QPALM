#include "lin_alg.h"
#include "global_opts.h"
#include <CUnit/CUnit.h>

#define DIM 3
#define TOL 1e-8

c_float a[DIM];
c_float b[DIM];
c_float c[DIM];

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
}

void test_vec_set_scalar_int(void){
    c_int scalar = 3;
    vec_set_scalar(a, scalar, DIM);
    for (size_t i = 0; i < DIM; i++) {
        CU_ASSERT_EQUAL(a[i], scalar);
    }
}

void test_vec_mult_scalar(void){
    c_int scalar = 3;
    vec_mult_scalar(a, scalar, DIM);
    CU_ASSERT_DOUBLE_EQUAL(a[0],0.3, TOL);
    CU_ASSERT_DOUBLE_EQUAL(a[1],7.5, TOL);
    CU_ASSERT_DOUBLE_EQUAL(a[2],-11.7, TOL);
}

void test_vec_prod(void){
    c_float expected[DIM+1] = {0.0,0.0,25,9.4};

    for (int i = 0; i == DIM; i++) {
        CU_ASSERT_DOUBLE_EQUAL(vec_prod(a, b, i), expected[i], TOL);
    }

}

void test_vec_add_scaled(void){
    int scalar = 4;
    vec_add_scaled(a, b, c, scalar, DIM);
    CU_ASSERT_DOUBLE_EQUAL(c[0], 0.1, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[1], 42.5, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[2], 12.1, TOL);
}

void test_vec_norm_inf(void){
    CU_ASSERT_DOUBLE_EQUAL(vec_norm_inf(a, DIM), 3.9, TOL);
    CU_ASSERT_DOUBLE_EQUAL(vec_norm_inf(b, DIM), 10, TOL);
}

void test_vec_ew_recipr(void){
    vec_ew_recipr(a, c, DIM);
    CU_ASSERT_DOUBLE_EQUAL(c[0], 10, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[1], 0.4, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[2], -0.256410256410256, TOL);
}

void test_vec_ew_max_vec(void){
    vec_ew_max_vec(a, b, c, DIM);
    CU_ASSERT_DOUBLE_EQUAL(c[0], 0.1, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[1], 10, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[2], 4, TOL);
}

void test_vec_ew_min_vec(void){
    vec_ew_min_vec(a, b, c, DIM);
    CU_ASSERT_DOUBLE_EQUAL(c[0], 0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[1], 2.5, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[2], -3.9, TOL);
}

void test_vec_ew_mid_vec(void){
    vec_ew_mid_vec(a, c, b, c, DIM);
    CU_ASSERT_DOUBLE_EQUAL(c[0], 0.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[1], 2.5, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[2], 0.0, TOL);
}

void test_vec_ew_prod(void){
    vec_ew_prod(a, b, c, DIM);
    CU_ASSERT_DOUBLE_EQUAL(c[0], 0.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[1], 25, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[2], -15.6, TOL);
}

void test_vec_ew_div(void){
    vec_ew_div(b, a, c, DIM);
    CU_ASSERT_DOUBLE_EQUAL(c[0], 0.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[1], 4, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[2], -1.025641025641026, TOL);
}

void test_vec_ew_sqrt(void){
    vec_ew_sqrt(b, c, DIM);
    CU_ASSERT_DOUBLE_EQUAL(c[0], 0.0, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[1], 3.162277660168380, TOL);
    CU_ASSERT_DOUBLE_EQUAL(c[2], 2.0, TOL);
}