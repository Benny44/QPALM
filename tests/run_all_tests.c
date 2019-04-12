#include <CUnit/Basic.h>
#include "test_lin_alg.h"

int main(){

    /* Initialize the CUnit test registry */
    if (CUE_SUCCESS != CU_initialize_registry())
        return CU_get_error();

    /* linear algebra tests */
    CU_TestInfo test_lin_alg[] = {
        { "test_vec_set_scalar", test_vec_set_scalar},
        { "test_vec_set_scalar_int", test_vec_set_scalar_int},
        CU_TEST_INFO_NULL,
    };

    /* list of suites to be tested */
    CU_SuiteInfo suites[] = {
        { "lin_alg", NULL, NULL, NULL, NULL, test_lin_alg},
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