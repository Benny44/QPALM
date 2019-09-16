#include "qpalm.h"
#include "global_opts.h"
#include "constants.h"
#include "load_data.h"

#define N 2000
#define M 3
#define ANZMAX 600
#define QNZMAX 380904

// static c_float solution[N] = {2.0000000e+00, -6.3801365e+01, -3.3821109e+03, -6.0483288e+00};


int main(void) {

    QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-4;
    settings->eps_rel = 1e-4;

    QPALMData *data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;

    cholmod_common common, *c;
    c = &common;
    CHOLMOD(start)(c);
    data->A = CHOLMOD(allocate_sparse)(data->m, data->n, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, c);
    data->Q = CHOLMOD(allocate_sparse)(data->n, data->n, QNZMAX , TRUE, TRUE, -1, CHOLMOD_REAL, c);
    CHOLMOD(finish)(c);
    
    load_data("../../../profiling/data/profiling_qp", data);

    // Setup workspace
    settings->verbose = FALSE;
    for(int k = 0; k < 1; k++) {
        QPALMWorkspace *work = qpalm_setup(data, settings, c);
        // Solve Problem

        qpalm_solve(work);

        c_print("Runtime: %.4f s\n", work->info->run_time);
        c_print("Iter: %ld\n", work->info->iter);
        //cleanup
        qpalm_cleanup(work);
    }

    c_free(settings);
    // Clean setup
    CHOLMOD(start)(c);
    CHOLMOD(free_sparse)(&data->Q, c);
    CHOLMOD(free_sparse)(&data->A, c);
    CHOLMOD(finish)(c);
    c_free(data->q);
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);
    return 0;
}