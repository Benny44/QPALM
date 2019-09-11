#include "minunit.h"
#include "test_lin_alg.h"
#include "test_basic_qp.h"
#include "test_prim_inf_qp.h"
#include "test_dua_inf_qp.h"
#include "test_degen_hess.h"
#include "test_cholmod_interface.h"
#include "test_nonconvex_qp.h"
#include "test_update.h"
#include "test_validate.h"
#include "test_error_handling.h"
#include "test_big_qp.h"
#include "test_ls_qp.h"


int main(){
    MU_INITIALIZE();
    MU_RUN_SUITE(suite_lin_alg);
    MU_RUN_SUITE(suite_cholmod);
    MU_RUN_SUITE(suite_basic_qp);
    MU_RUN_SUITE(suite_prim_inf_qp);
    MU_RUN_SUITE(suite_dua_inf_qp);
    MU_RUN_SUITE(suite_degen_hess);
    MU_RUN_SUITE(suite_nonconvex);
    MU_RUN_SUITE(suite_update);
    MU_RUN_SUITE(suite_validation);
    MU_RUN_SUITE(suite_error_handling);
    MU_RUN_SUITE(suite_big_qp);
    MU_RUN_SUITE(suite_ls_qp);
    MU_REPORT();
    
    return minunit_fail; /* =0 if all tests passed, >0 otherwise */
	    
}