#include "lin_alg.h"
#include "global_opts.h"
#include "minunit.h"

#define DIM 3
#define TOL 1e-8

c_float a[DIM];
c_int a_int[DIM];
c_float b[DIM];
c_float c[DIM];

void reset_abc(void) {
    a[0] = 0.1; a[1] = 2.5; a[2] = -3.9;
    b[0] = 0.0; b[1] = 10; b[2] = 4;
    c[0] = 0.0; c[1] = 0.0; c[2] = 0.0;
}

MU_TEST(test_vec_set_scalar) {
    c_float scalar = 5.5;
    vec_set_scalar(a, scalar, DIM);
    for (size_t i = 0; i < DIM; i++) {
        mu_assert_double_eq(a[i], scalar, TOL);
    }
}

MU_TEST(test_vec_set_scalar_int){
    c_int scalar = 3;
    vec_set_scalar_int(a_int, scalar, DIM);
    for (size_t i = 0; i < DIM; i++) {
        mu_assert_long_eq(a_int[i], scalar);
    }
}

MU_TEST(test_vec_self_mult_scalar){
    c_float scalar = 3.0;
    vec_self_mult_scalar(a, scalar, DIM);
    mu_assert_double_eq(a[0], 0.3, TOL);
    mu_assert_double_eq(a[1], 7.5, TOL);
    mu_assert_double_eq(a[2], -11.7, TOL);
}

MU_TEST(test_vec_prod){
    c_float expected[DIM+1] = {0.0,0.0,25,9.4};

    for (size_t i = 0; i == DIM; i++) {
        mu_assert_double_eq(vec_prod(a, b, i), expected[i], TOL);
    }

}

MU_TEST(test_vec_add_scaled){
    int scalar = 4;
    vec_add_scaled(a, b, c, scalar, DIM);
    mu_assert_double_eq(c[0], 0.1, TOL);
    mu_assert_double_eq(c[1], 42.5, TOL);
    mu_assert_double_eq(c[2], 12.1, TOL);
}

MU_TEST(test_vec_norm_inf){
    mu_assert_double_eq(vec_norm_inf(a, DIM), 3.9, TOL);
    mu_assert_double_eq(vec_norm_inf(b, DIM), 10, TOL);
}

MU_TEST(test_vec_ew_recipr){
    vec_ew_recipr(a, c, DIM);
    mu_assert_double_eq(c[0], 10, TOL);
    mu_assert_double_eq(c[1], 0.4, TOL);
    mu_assert_double_eq(c[2], -0.256410256410256, TOL);
}

MU_TEST(test_vec_ew_max_vec){
    vec_ew_max_vec(a, b, c, DIM);
    mu_assert_double_eq(c[0], 0.1, TOL);
    mu_assert_double_eq(c[1], 10, TOL);
    mu_assert_double_eq(c[2], 4, TOL);
}

MU_TEST(test_vec_ew_min_vec){
    vec_ew_min_vec(a, b, c, DIM);
    mu_assert_double_eq(c[0], 0, TOL);
    mu_assert_double_eq(c[1], 2.5, TOL);
    mu_assert_double_eq(c[2], -3.9, TOL);
}

MU_TEST(test_vec_ew_mid_vec){
    vec_ew_mid_vec(a, c, b, c, DIM);
    mu_assert_double_eq(c[0], 0.0, TOL);
    mu_assert_double_eq(c[1], 2.5, TOL);
    mu_assert_double_eq(c[2], 0.0, TOL);
}

MU_TEST(test_vec_ew_prod){
    vec_ew_prod(a, b, c, DIM);
    mu_assert_double_eq(c[0], 0.0, TOL);
    mu_assert_double_eq(c[1], 25, TOL);
    mu_assert_double_eq(c[2], -15.6, TOL);
}

MU_TEST(test_vec_ew_div){
    vec_ew_div(b, a, c, DIM);
    mu_assert_double_eq(c[0], 0.0, TOL);
    mu_assert_double_eq(c[1], 4, TOL);
    mu_assert_double_eq(c[2], -1.025641025641026, TOL);
}

MU_TEST(test_vec_ew_sqrt){
    vec_ew_sqrt(b, c, DIM);
    mu_assert_double_eq(c[0], 0.0, TOL);
    mu_assert_double_eq(c[1], 3.162277660168380, TOL);
    mu_assert_double_eq(c[2], 2.0, TOL);
}

MU_TEST_SUITE(suite_lin_alg) {
    MU_SUITE_CONFIGURE(NULL, NULL, reset_abc, NULL);

    MU_RUN_TEST(test_vec_set_scalar);
    MU_RUN_TEST(test_vec_set_scalar_int);
    MU_RUN_TEST(test_vec_self_mult_scalar);
    MU_RUN_TEST(test_vec_prod);
    MU_RUN_TEST(test_vec_add_scaled);
    MU_RUN_TEST(test_vec_norm_inf);
    MU_RUN_TEST(test_vec_ew_recipr);
    MU_RUN_TEST(test_vec_ew_max_vec);
    MU_RUN_TEST(test_vec_ew_min_vec);
    MU_RUN_TEST(test_vec_ew_mid_vec);
    MU_RUN_TEST(test_vec_ew_prod);
    MU_RUN_TEST(test_vec_ew_div);
    MU_RUN_TEST(test_vec_ew_sqrt);
    
}