#include "types.h"
#include "util.h"
#include "constants.h"
#include "global_opts.h"
#include "load_data.h"
#include "cholmod_interface.h"

#define N 2000
#define M 3
#define ANZMAX 600
// #define QNZMAX 380904
#define QNZMAX 2000

#   include <time.h>
#   include <sys/time.h>
# include <stdio.h>
#include <stdlib.h>

int main(void) {

    QPALMData *data = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = N;
    data->m = M;

    cholmod_common common, *c;
    c = &common;
    CHOLMOD(start)(c);
    data->A = CHOLMOD(allocate_sparse)(data->m, data->n, ANZMAX, TRUE, TRUE, 0, CHOLMOD_REAL, c);
    data->Q = CHOLMOD(allocate_sparse)(data->n, data->n, QNZMAX , TRUE, TRUE, -1, CHOLMOD_REAL, c);

    load_data("../../../profiling/data/profiling_qp", data);

    QPALMTimer* timer = c_malloc(sizeof(QPALMTimer));

    c_float runtime;
    
    cholmod_dense *y = CHOLMOD(allocate_dense)(M, 1, M, CHOLMOD_REAL, c);
    cholmod_dense *Aty = CHOLMOD(allocate_dense)(N, 1, N, CHOLMOD_REAL, c);
    cholmod_dense *x = CHOLMOD(allocate_dense)(N, 1, N, CHOLMOD_REAL, c);

    c_float *yx = y->x;
    
    for(c_int i =0; i < M; i++) {
        yx[i] = (c_float) rand()/RAND_MAX;
    }


    qpalm_tic(timer);
    for (c_int i = 0; i < 1000; i++) {
        mat_tpose_vec(data->A, y, Aty, c);
    }
    runtime = qpalm_toc(timer);   

    printf("At*y runtime: %.4e \n", runtime);

    qpalm_tic(timer);
    for (c_int i = 0; i < 100; i++) {
        mat_tpose_vec(data->Q, Aty, x, c);
    }
    runtime = qpalm_toc(timer);   

    printf("Q*x runtime: %.4e \n", runtime);
    
    // Clean setup
    CHOLMOD(free_sparse)(&data->Q, c);
    CHOLMOD(free_sparse)(&data->A, c);
    CHOLMOD(finish)(c);
    c_free(data->q);
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);
    return 0;
}