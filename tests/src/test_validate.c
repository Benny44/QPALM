#include "minunit.h"
#include "validate.h"
#include "global_opts.h"
#include "constants.h"
#include "types.h"
#include "test_validate.h"
#include "qpalm.h"
#include <stdio.h>

#define M 3

QPALMSettings *settings;
QPALMData *data;

void validate_suite_setup(void) {

    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->m = M;

    data->bmin = c_calloc(M,sizeof(c_float));
    data->bmin[0] = -1.0; data->bmin[1] = -1.0; data->bmin[2] = -1.0; 
    data->bmax = c_calloc(M,sizeof(c_float));
    data->bmax[0] = 1.0; data->bmax[1] = 1.0; data->bmax[2] = 1.0; 
}

void validate_suite_teardown(void) {
    c_free(settings);
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);
}

void validate_test_setup(void) {
    if (data==NULL) {
        data = (QPALMData *)c_malloc(sizeof(QPALMData));
        data->m = M;
        data->bmin = c_calloc(M,sizeof(c_float));
        data->bmax = c_calloc(M,sizeof(c_float));
    }
    if (settings == NULL) {
        settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    }
    
    qpalm_set_default_settings(settings);

    data->bmin[0] = -1.0; data->bmin[1] = -1.0; data->bmin[2] = -1.0;
    data->bmax[0] = 1.0; data->bmax[1] = 1.0; data->bmax[2] = 1.0;
}



/* Data validation */
MU_TEST(test_correct_data) {
    mu_assert_true(validate_data(data));
}

MU_TEST(test_missing_data) {
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);
    data=NULL;    
    mu_assert_false(validate_data(data));    
}

MU_TEST(test_bounds_mismatch) {
    data->bmin[1] = 2; //higher than the upper bound
    mu_assert_false(validate_data(data));
}

/*Settings validation */
MU_TEST(test_correct_settings) {
    mu_assert_true(validate_settings(settings));
}

MU_TEST(test_missing_settings) {
    c_free(settings);
    settings=NULL;
    mu_assert_false(validate_settings(settings));
}
MU_TEST(test_max_iter_out_of_bounds){
    settings->max_iter = -1;
    mu_assert_false(validate_settings(settings));
    settings->max_iter = 0;
    mu_assert_false(validate_settings(settings));
    settings->max_iter = 1;
    mu_assert_true(validate_settings(settings));
}
MU_TEST(test_tol_out_of_bounds){
    settings->eps_abs = -1;
    mu_assert_false(validate_settings(settings));
    settings->eps_abs = 1;

    settings->eps_rel = -1;
    mu_assert_false(validate_settings(settings));

    settings->eps_rel = 0; settings->eps_abs = 0;
    mu_assert_false(validate_settings(settings));
    settings->eps_abs = 1;

    settings->eps_abs_in = -1;
    mu_assert_false(validate_settings(settings));
    settings->eps_abs_in = 1;

    settings->eps_rel_in = -1;
    mu_assert_false(validate_settings(settings));
    settings->eps_rel_in = 0; settings->eps_abs_in = 0;
    mu_assert_false(validate_settings(settings));
    settings->eps_abs_in = 1;

    settings->eps_dual_inf = -1;
    mu_assert_false(validate_settings(settings));
    settings->eps_dual_inf = 1;

    settings->eps_prim_inf = -1;
    mu_assert_false(validate_settings(settings));
}

MU_TEST(test_rho_out_of_bounds){
    settings->rho = 0;
    mu_assert_false(validate_settings(settings));
    settings->rho = 1;
    mu_assert_false(validate_settings(settings));
    settings->rho = 0.5;
    mu_assert_true(validate_settings(settings));
}

MU_TEST(test_theta_out_of_bounds){
    settings->theta = 1;
    mu_assert_true(validate_settings(settings));
    settings->theta = 2;
    mu_assert_false(validate_settings(settings));
}

MU_TEST(test_delta_out_of_bounds){
    settings->delta = 0.5;
    mu_assert_false(validate_settings(settings));
    settings->delta = 1;
    mu_assert_false(validate_settings(settings));
    settings->delta = 1.1;
    mu_assert_true(validate_settings(settings));
}

MU_TEST(test_gamma_out_of_bounds){
    settings->gamma_init = 0;
    mu_assert_false(validate_settings(settings));
    settings->gamma_init = 1e5;
    
    settings->gamma_upd = 0.5;
    mu_assert_false(validate_settings(settings));
    settings->gamma_upd = 1;
    mu_assert_true(validate_settings(settings));

    settings->gamma_max = settings->gamma_init;
    mu_assert_true(validate_settings(settings));
    settings->gamma_max /= 10;
    mu_assert_false(validate_settings(settings));
}

MU_TEST(test_scaling_out_of_bounds){
    settings->scaling = -1;
    mu_assert_false(validate_settings(settings));
    settings->scaling = 0;
    mu_assert_true(validate_settings(settings));
}

MU_TEST(test_booleans){
    settings->proximal = FALSE;
    mu_assert_true(validate_settings(settings));
    settings->proximal = 3;
    mu_assert_false(validate_settings(settings));
    settings->proximal = TRUE;

    settings->warm_start = FALSE;
    mu_assert_true(validate_settings(settings));
    settings->warm_start = 3;
    mu_assert_false(validate_settings(settings));
    settings->warm_start = TRUE;

    settings->verbose = FALSE;
    mu_assert_true(validate_settings(settings));
    settings->verbose = 3;
    mu_assert_false(validate_settings(settings));
    settings->verbose = TRUE;
}

MU_TEST_SUITE(suite_validation) {
    MU_SUITE_CONFIGURE(validate_suite_setup, validate_suite_teardown, validate_test_setup, NULL);

    MU_RUN_TEST(test_correct_data);
    MU_RUN_TEST(test_missing_data);
    MU_RUN_TEST(test_bounds_mismatch);
    MU_RUN_TEST(test_correct_settings);
    MU_RUN_TEST(test_missing_settings);
    MU_RUN_TEST(test_max_iter_out_of_bounds);
    MU_RUN_TEST(test_tol_out_of_bounds);
    MU_RUN_TEST(test_rho_out_of_bounds);
    MU_RUN_TEST(test_theta_out_of_bounds);
    MU_RUN_TEST(test_delta_out_of_bounds);
    MU_RUN_TEST(test_gamma_out_of_bounds);
    MU_RUN_TEST(test_scaling_out_of_bounds);
    MU_RUN_TEST(test_booleans);
}