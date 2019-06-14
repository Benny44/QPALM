#include "validate.h"
#include "global_opts.h"
#include "constants.h"
#include "types.h"
#include "test_validate.h"
#include "qpalm.h"
#include <stdio.h>

#include <CUnit/CUnit.h>

#define M 3

QPALMSettings *settings;
QPALMData *data;

int validate_suite_setup(void) {

    settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->m = M;
    // c_float bmin[3] = {-1.0, -1.0, -1.0};
    // c_float bmax[3] = {1.0, 1.0, 1.0};
    // data->bmin = bmin;
    // data->bmax = bmax;

    data->bmin = c_calloc(M,sizeof(c_float));
    data->bmin[0] = -1.0; data->bmin[1] = -1.0; data->bmin[2] = -1.0; 
    data->bmax = c_calloc(M,sizeof(c_float));
    data->bmax[0] = 1.0; data->bmax[1] = 1.0; data->bmax[2] = 1.0; 


    return 0;
}

int validate_suite_teardown(void) {
    c_free(settings);
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);

    return 0;
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
void test_correct_data(void) {
    CU_ASSERT_TRUE(validate_data(data));
}

void test_missing_data(void) {
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);
    data=NULL;    
    CU_ASSERT_FALSE(validate_data(data));    
}

void test_bounds_mismatch(void) {
    data->bmin[1] = 2; //higher than the upper bound
    CU_ASSERT_FALSE(validate_data(data));
}

/*Settings validation */
void test_correct_settings(void) {
    CU_ASSERT_TRUE(validate_settings(settings));
}

void test_missing_settings(void) {
    c_free(settings);
    settings=NULL;
    CU_ASSERT_FALSE(validate_settings(settings));
}
void test_max_iter_out_of_bounds(void){
    settings->max_iter = -1;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->max_iter = 0;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->max_iter = 1;
    CU_ASSERT_TRUE(validate_settings(settings));
}
void test_tol_out_of_bounds(void){
    settings->eps_abs = -1;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->eps_abs = 1;

    settings->eps_rel = -1;
    CU_ASSERT_FALSE(validate_settings(settings));

    settings->eps_rel = 0; settings->eps_abs = 0;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->eps_abs = 1;

    settings->eps_abs_in = -1;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->eps_abs_in = 1;

    settings->eps_rel_in = -1;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->eps_rel_in = 0; settings->eps_abs_in = 0;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->eps_abs_in = 1;

    settings->eps_dual_inf = -1;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->eps_dual_inf = 1;

    settings->eps_prim_inf = -1;
    CU_ASSERT_FALSE(validate_settings(settings));
}

void test_rho_out_of_bounds(void){
    settings->rho = 0;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->rho = 1;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->rho = 0.5;
    CU_ASSERT_TRUE(validate_settings(settings));
}

void test_theta_out_of_bounds(void){
    settings->theta = 1;
    CU_ASSERT_TRUE(validate_settings(settings));
    settings->theta = 2;
    CU_ASSERT_FALSE(validate_settings(settings));
}

void test_delta_out_of_bounds(void){
    settings->delta = 0.5;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->delta = 1;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->delta = 1.1;
    CU_ASSERT_TRUE(validate_settings(settings));
}

void test_gamma_out_of_bounds(void){
    settings->gamma_init = 0;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->gamma_init = 1e5;
    
    settings->gamma_upd = 0.5;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->gamma_upd = 1;
    CU_ASSERT_TRUE(validate_settings(settings));

    settings->gamma_max = settings->gamma_init;
    CU_ASSERT_TRUE(validate_settings(settings));
    settings->gamma_max /= 10;
    CU_ASSERT_FALSE(validate_settings(settings));
}

void test_scaling_out_of_bounds(void){
    settings->scaling = -1;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->scaling = 0;
    CU_ASSERT_TRUE(validate_settings(settings));
}

void test_booleans(void){
    settings->proximal = FALSE;
    CU_ASSERT_TRUE(validate_settings(settings));
    settings->proximal = 3;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->proximal = TRUE;

    settings->warm_start = FALSE;
    CU_ASSERT_TRUE(validate_settings(settings));
    settings->warm_start = 3;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->warm_start = TRUE;

    settings->verbose = FALSE;
    CU_ASSERT_TRUE(validate_settings(settings));
    settings->verbose = 3;
    CU_ASSERT_FALSE(validate_settings(settings));
    settings->verbose = TRUE;
}
