#include "qpalm.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#ifdef USE_LADEL
#include "ladel.h"
#elif defined USE_CHOLMOD
#include "cholmod.h"
#endif

solver_sparse* mtx_load_A(FILE *fp, size_t *n, size_t *m){
    // Get the first line out of the way
    char first_line[1000];
    if (fgets(first_line, sizeof first_line, fp) == NULL) {
        fprintf(stderr, "Empty file\n");
    };

    // Read the matrix size
    size_t nnz;
    if (fscanf(fp, "%lu %lu %lu", m, n, &nnz) != 3) {
        fprintf(stderr, "Wrong file format. Expected second line to contain row col nnz.\n");
    };
    
    solver_common c;
    solver_sparse *A;
    #ifdef USE_LADEL
    A = ladel_sparse_alloc(*m, *n, nnz, UNSYMMETRIC, TRUE, FALSE);
    #elif defined USE_CHOLMOD
    CHOLMOD(start)(&c);
    A = CHOLMOD(allocate_sparse)(*m, *n, nnz, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    CHOLMOD(finish)(&c);
    #endif

    c_float *Ax; c_int *Ai, *Ap;
    Ax = A->x; Ap = A->p; Ai = A->i;
    Ap[0] = 0;

    size_t elem = 0, col = 1, prev_col = 1;
    c_int row;
    c_float temp;

    while (fscanf(fp, "%lu %lu %le", &row, &col, &temp) != EOF) {

        Ai[elem] = --row;
        
        if(col > prev_col) {
            for (; prev_col < col; prev_col++) {
                Ap[prev_col+1] = Ap[prev_col];
            }
        }
        Ap[col]++;

        if (temp > QPALM_INFTY) {
            temp = QPALM_INFTY;
        } else if (temp < -QPALM_INFTY) {
            temp = -QPALM_INFTY;
        }
        Ax[elem] = temp;

        elem++;
    }
    

    return A;
}

solver_sparse* mtx_load_Q(FILE *fp, size_t n_check){
        // Get the first line out of the way
    char first_line[1000];
    if (fgets(first_line, sizeof first_line, fp) == NULL){
        fprintf(stderr, "Empty file\n");
    }

    // Read the matrix size
    size_t m, n, nnz;
    if (fscanf(fp, "%lu %lu %lu", &m, &n, &nnz) != 3) {
        fprintf(stderr, "Wrong file format. Expected second line to contain row col nnz.\n");
    };

    if ((m != n) || (n != n_check)) {
        fprintf(stderr, "Expected a square matrix of %lu by %lu but got a matrix of %lu by %lu\n", n_check, n_check, m, n);
    }
    
    solver_common c;
    solver_sparse *Q;
    #ifdef USE_LADEL
    Q = ladel_sparse_alloc(n, n, nnz, UPPER, TRUE, FALSE);
    #elif defined USE_CHOLMOD
    CHOLMOD(start)(&c);
    Q = CHOLMOD(allocate_sparse)(n, n, nnz, TRUE, TRUE, -1, CHOLMOD_REAL, &c);
    CHOLMOD(finish)(&c);
    #endif

    c_float *Qx; c_int *Qi, *Qp;
    Qx = Q->x; Qp = Q->p; Qi = Q->i;
    Qp[0] = 0;

    size_t elem = 0, col = 1, prev_col = 1;
    c_int row;
    c_float temp;

    while (fscanf(fp, "%lu %lu %le", &row, &col, &temp) != EOF) {

        Qi[elem] = --row;
        
        if (col > prev_col) {
            for (; prev_col < col; prev_col++) {
                Qp[prev_col+1] = Qp[prev_col];
            }          
        }
        Qp[col]++;
        if (temp > QPALM_INFTY) {
            temp = QPALM_INFTY;
        } else if (temp < -QPALM_INFTY) {
            temp = -QPALM_INFTY;
        }
        Qx[elem] = temp;

        elem++;
    }

    return Q;
}

c_float* mtx_load_dense(FILE* fp, size_t n_check) {
    // Get the first line out of the way
    char first_line[1000];
    if (fgets(first_line, sizeof first_line, fp) == NULL){
        fprintf(stderr, "Empty file\n");
    }

    // Read the vector sizes
    size_t n, col, nnz;
    if (fscanf(fp, "%lu %lu %lu", &n, &col, &nnz) != 3) {
        fprintf(stderr, "Wrong file format. Expected second line to contain row col nnz.\n");
    };
    if (col != 1 || n != n_check) {
        fprintf(stderr, "Expected a column vector of length %lu but got a matrix of %lu by %lu\n", n_check, n, col);
    }
    c_float *result = (c_float *) calloc(n, sizeof(c_float));
    size_t i;
    c_float temp;
    while (fscanf(fp, "%lu %*u %le", &i, &temp) != EOF) {
        if (temp > QPALM_INFTY) {
            temp = QPALM_INFTY;
        } else if (temp < -QPALM_INFTY) {
            temp = -QPALM_INFTY;
        }
        result[--i] = temp;
    }
    return result;

}

int main(int argc, char*argv[]){
    
    if (argc != 6) {
        fprintf(stderr, "Wrong number of arguments. Correct usage is qpalm_mtx A.mtx Q.mtx q.mtx lba.mtx uba.mtx\n");
    }

    // Load problem data
    size_t n, m, fileno;
    FILE* fp;

    fileno = 0;
    // Load A (plus determine n and m from A)
    fileno++;
    fp = fopen(argv[fileno], "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", argv[fileno]);return 1;
    }
    solver_sparse* A = mtx_load_A(fp, &n, &m);
    fclose(fp);

    // Load Q
    fileno++;
    fp = fopen(argv[fileno], "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", argv[fileno]);return 1;
    }
    solver_sparse* Q = mtx_load_Q(fp, n);
    fclose(fp);

    // Load q
    fileno++;
    fp = fopen(argv[fileno], "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", argv[fileno]);return 1;
    }
    c_float* q = mtx_load_dense(fp, n);
    fclose(fp);

    // Load lba
    fileno++;
    fp = fopen(argv[fileno], "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", argv[fileno]);return 1;
    }
    c_float* bmin = mtx_load_dense(fp, m);
    fclose(fp);

    // Load uba
    fileno++;
    fp = fopen(argv[fileno], "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", argv[fileno]);return 1;
    }
    c_float* bmax = mtx_load_dense(fp, m);
    fclose(fp);

    // Problem settings
    QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));

    // Structures
    QPALMWorkspace *work; // Workspace
    QPALMData *data;      // QPALMData

    // Populate data
    data    = (QPALMData *)c_malloc(sizeof(QPALMData));
    data->n = n;
    data->m = m;
    data->c = 0;
    data->q = q;
    data->bmin = bmin;
    data->bmax = bmax;
    data->A = A;
    data->Q = Q;
    

    // Define Solver settings as default
    qpalm_set_default_settings(settings);

    // Setup workspace
    solver_common c;
    work = qpalm_setup(data, settings);

    // Solve Problem
    qpalm_solve(work);

    // Clean workspace
    #ifdef USE_LADEL
    data->A = ladel_sparse_free(data->A);
    data->Q = ladel_sparse_free(data->Q);
    #elif defined USE_CHOLMOD
    CHOLMOD(start)(&c);
    CHOLMOD(free_sparse)(&data->Q, &c);
    CHOLMOD(free_sparse)(&data->A, &c);
    CHOLMOD(finish)(&c);
    #endif

    qpalm_cleanup(work);
    c_free(data->q);
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);
    c_free(settings);

    return 0;
}