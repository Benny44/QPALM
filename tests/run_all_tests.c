#include <CUnit/Basic.h>
#include "test_lin_alg.h"
#include "test_basic_qp.h"
#include "test_prim_inf_qp.h"
#include "test_dua_inf_qp.h"
#include "test_cholmod_interface.h"
#include "test_validate.h"

int main(){

    /* Initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* linear algebra tests */
    CU_TestInfo suite_lin_alg[] = {
        { "test_vec_set_scalar", test_vec_set_scalar},
        { "test_vec_set_scalar_int", test_vec_set_scalar_int},
        { "test_vec_mult_scalar", test_vec_mult_scalar},
        { "test_vec_prod", test_vec_prod},
        { "test_vec_add_scaled", test_vec_add_scaled},
        { "test_vec_norm_inf", test_vec_norm_inf},
        { "test_vec_ew_recipr", test_vec_ew_recipr},
        { "test_vec_ew_max_vec", test_vec_ew_max_vec},
        { "test_vec_ew_min_vec", test_vec_ew_min_vec},
        { "test_vec_ew_mid_vec", test_vec_ew_mid_vec},
        { "test_vec_ew_prod", test_vec_ew_prod},
        { "test_vec_ew_div", test_vec_ew_div},
        { "test_vec_ew_sqrt", test_vec_ew_sqrt},
        CU_TEST_INFO_NULL,
    };

    /* cholmod interface tests */
    CU_TestInfo suite_cholmod[] = {
        { "test_mat_vec", test_mat_vec},
        { "test_mat_tpose_vec", test_mat_tpose_vec},
        { "test_mat_inf_norm_cols", test_mat_inf_norm_cols},
        { "test_mat_inf_norm_rows", test_mat_inf_norm_rows},
        { "test_ldlchol", test_ldlchol},
        CU_TEST_INFO_NULL,
    };

    /* basic qp */
    CU_TestInfo suite_basic_qp[] = {
        { "test_basic_qp", test_basic_qp},
        { "test_basic_qp_unscaled", test_basic_qp_unscaled},
        { "test_basic_qp_noprox", test_basic_qp_noprox},
        { "test_basic_qp_noprox_unscaled", test_basic_qp_noprox_unscaled},
        { "test_basic_qp_warm_start", test_basic_qp_warm_start},     
        { "test_basic_qp_warm_start_unscaled", test_basic_qp_warm_start_unscaled},     
        { "test_basic_qp_warm_start_noprox", test_basic_qp_warm_start_noprox},     
        { "test_basic_qp_warm_start_noprox_unscaled", test_basic_qp_warm_start_noprox_unscaled},     
        CU_TEST_INFO_NULL,
    };

    /* primal infeasible qp */
    CU_TestInfo suite_prim_inf_qp[] = {
        { "test_prim_inf_qp", test_prim_inf_qp},
        { "test_prim_inf_qp_unscaled", test_prim_inf_qp_unscaled},
        { "test_prim_inf_qp_noprox", test_prim_inf_qp_noprox},
        { "test_prim_inf_qp_noprox_unscaled", test_prim_inf_qp_noprox_unscaled},
        CU_TEST_INFO_NULL,
    };
    
    /* primal infeasible qp */
    CU_TestInfo suite_dua_inf_qp[] = {
        { "test_dua_inf_qp", test_dua_inf_qp},
        { "test_dua_inf_qp_unscaled", test_dua_inf_qp_unscaled},
        { "test_dua_inf_qp_noprox", test_dua_inf_qp_noprox},
        { "test_dua_inf_qp_noprox_unscaled", test_dua_inf_qp_noprox_unscaled},
        CU_TEST_INFO_NULL,
    };

    /* validation functions */
    CU_TestInfo suite_validation[] = {
        { "test_correct_data", test_correct_data},
        { "test_missing_data", test_missing_data},
        { "test_bounds_mismatch", test_bounds_mismatch},
        { "test_correct_settings", test_correct_settings},
        { "test_missing_settings", test_missing_settings},
        { "test_max_iter_out_of_bounds", test_max_iter_out_of_bounds},
        { "test_tol_out_of_bounds", test_tol_out_of_bounds},
        { "test_rho_out_of_bounds", test_rho_out_of_bounds},
        { "test_theta_out_of_bounds", test_theta_out_of_bounds},
        { "test_delta_out_of_bounds", test_delta_out_of_bounds},
        { "test_gamma_out_of_bounds", test_gamma_out_of_bounds},
        { "test_scaling_out_of_bounds", test_scaling_out_of_bounds},
        { "test_booleans", test_booleans},        
        CU_TEST_INFO_NULL,
    };
    
    /* list of suites to be tested */
    CU_SuiteInfo suites[] = {
        { "lin_alg", NULL, NULL, reset_abc, NULL, suite_lin_alg},
        { "cholmod", cholmod_qp_setup, cholmod_qp_teardown, cholmod_set_QdAd, NULL, suite_cholmod},
        { "basic_qp", basic_qp_suite_setup, basic_qp_suite_teardown, NULL, basic_qp_test_teardown, suite_basic_qp},
        { "prim_inf_qp", prim_inf_qp_suite_setup, prim_inf_qp_suite_teardown, NULL, prim_inf_qp_test_teardown, suite_prim_inf_qp},
        { "dua_inf_qp", dua_inf_qp_suite_setup, dua_inf_qp_suite_teardown, NULL, dua_inf_qp_test_teardown,suite_dua_inf_qp},
        { "validation", validate_suite_setup, validate_suite_teardown, validate_test_setup, NULL, suite_validation},
        CU_SUITE_INFO_NULL,
    };

    if (CUE_SUCCESS != CU_register_suites(suites)) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* Run all tests using the CUnit Basic interface */
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    CU_cleanup_registry();
    return CU_get_error();
        

}