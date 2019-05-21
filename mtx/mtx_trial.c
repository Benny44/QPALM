#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>


typedef double c_float;

c_float* mtx_load_dense(FILE* fp) {
    // Get the first line out of the way
    char first_line[1000];
    fgets(first_line, sizeof first_line, fp);

    // Read the vector sizes
    size_t n, col, nnz;
    fscanf(fp, "%lu %lu %lu", &n, &col, &nnz);
    printf("n: %lu\n",n);
    if (col != 1) {
        fprintf(stderr, "Expected a column vector but got a matrix of %ld by %ld\n", n, col);
    }
    c_float *result = (c_float *) calloc(n, sizeof(c_float));
    result[10] = 4;
    size_t i;
    c_float temp;
    int fscanfres;
    while ((fscanfres = fscanf(fp, "%lu %*u %le", &i, &temp)) != EOF) {
        if (temp > 1e20) {
            temp = 1e20;
        } else if (temp < -1e20) {
            temp = -1e20;
        }
        result[--i] = temp;
    }
    return result;
}

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Wrong number of arguments. Correct usage is \n");
    }

    size_t n, m, fileno;
    FILE* fp;

    fileno = 0;
    // Load q
    fileno++;
    fp = fopen(argv[fileno], "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", argv[fileno]);
    }

    printf("Before loading dense\n");
    c_float* q = mtx_load_dense(fp);

    printf("resulting q: ");
    for (int i = 0; i < 12; i++) {
        printf(" %f", q[i]);
    }
    printf("\n");
}