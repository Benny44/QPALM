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
    MU_REPORT();
    //if (minunit_fail) {
    //    return 1;
    //} else {
    return minunit_fail;
    //}
	    
}