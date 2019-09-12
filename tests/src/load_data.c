#include "types.h"
#include "global_opts.h"
#include "constants.h"
#include "cholmod.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

void load_sparse_matrix(FILE *fp, cholmod_sparse *A){

    c_int m = A->nrow;
    c_int n = A->ncol;

    c_float *Ax; c_int *Ai, *Ap;
    Ax = A->x; Ap = A->p; Ai = A->i;
    Ap[0] = 0; 

    c_float temp;
    c_int row, col, prev_col, k, elem; k = elem = 0; prev_col = 0;

    while (fscanf(fp, "%le", &temp) != EOF) {
        
        if (temp != 0) {
            row = mod(k,m); col = (k/m)+1;
            Ai[elem] = row;
            Ax[elem] = temp;

            if (col > prev_col) {
                for (; prev_col < col; prev_col++) {
                    Ap[prev_col+1] = Ap[prev_col];
                }          
            }
            Ap[col]++;   
            elem++;
        }
        k++;
        
    }

    for (; col < n; col++) {
        Ap[col+1] = Ap[col];
    }

    c_print("in A/Q: nnz = %ld, size = %ld\n", elem, k);
    for (k = 0; k < elem; k++) {
        c_print("Ax[%ld] = %.16f;\n", k, Ax[k]);
        c_print("Ai[%ld] = %ld;\n", k, Ai[k]);
    }
    for (k = 0; k < n+1; k++) {
        c_print("Ap[%ld] = %ld;\n", k, Ap[k]);
    }
}

c_float *load_dense(FILE *fp, size_t n) {
    c_float temp;
    c_float *result = c_calloc(n, sizeof(c_float));
    size_t k = 0;
     while (fscanf(fp, "%le", &temp) != EOF) {
         result[k] = temp;
         k++;
     }
     return result;
}

void load_data(const char* name, QPALMData* data) {
    

    char* name_A = c_calloc(strlen(name)+2+1, sizeof(char)); 
    name_A[0] = '\0';
    strcat(name_A, name);
    strcat(name_A, "_A");

    char* name_Q = c_calloc(strlen(name)+2+1, sizeof(char));
    name_Q[0] = '\0';
    strcat(name_Q, name);
    strcat(name_Q, "_Q");

    char* name_q = c_calloc(strlen(name)+2+1, sizeof(char));
    name_q[0] = '\0';
    strcat(name_q, name);
    strcat(name_q, "_q");

    char* name_bmin = c_calloc(strlen(name)+5+1, sizeof(char));
    name_bmin[0] = '\0';
    strcat(name_bmin, name);
    strcat(name_bmin, "_bmin");

    char* name_bmax = c_calloc(strlen(name)+5+1, sizeof(char));
    name_bmax[0] = '\0';
    strcat(name_bmax, name);
    strcat(name_bmax, "_bmax");


    FILE *fp;
    fp = fopen(name_A, "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", name_A);
    }
    load_sparse_matrix(fp, data->A);
    fclose(fp);

    fp = fopen(name_Q, "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", name_Q);
    }
    load_sparse_matrix(fp, data->Q);
    fclose(fp);

    fp = fopen(name_q, "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", name_q);
    }
    data->q = load_dense(fp, data->n);
    fclose(fp);

    fp = fopen(name_bmin, "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", name_bmin);
    }
    data->bmin = load_dense(fp, data->m);
    fclose(fp);

    fp = fopen(name_bmax, "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", name_bmax);
    }
    data->bmax = load_dense(fp, data->m);
    fclose(fp);

    c_free(name_A);
    c_free(name_Q);
    c_free(name_q);
    c_free(name_bmin);
    c_free(name_bmax);

}