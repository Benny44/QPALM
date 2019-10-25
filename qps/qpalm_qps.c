#include "qpalm.h"
#include "constants.h"
#include "global_opts.h"
#include "cholmod.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char*argv[]){

    if (argc != 2) {
        fprintf(stderr, "Wrong number of arguments. Correct usage is qpalm_qps problem.qps\n");
    }
    
    // Load problem data
    size_t n, m;
    FILE* fp;

    fp = fopen(argv[1], "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", argv[1]);
        return 1;
    }

    char line[100], command[20], name[50];
    fgets(line, 100, fp);
    if (sscanf(line, "%s %s", command, name) != 2) {
        fprintf(stderr, "Wrong file format. Expected first line to contain NAME problem_name nnz.\n");
        return 1;
    };

    if (strcmp(command, "NAME")) {
        fprintf(stderr, "Wrong file format. Expected first line to contain NAME problem_name nnz.\n");
        return 1;
    }
    printf("Reading problem %s\n", name);

    fgets(line, 100, fp);
    sscanf(line, "%s", command);
    if (strcmp(command, "ROWS")) {
        fprintf(stderr, "Wrong file format. Expected next command to be ROWS.\n");
        return 1;
    }

    char next_char;
    next_char = fgetc(fp);
    char NLGE[1], buf[20], objective[20];
    m=0; //keep track of number of constraints
    char constraint_signs[100000];

    while(next_char == ' ') {
        fgets(line, 100, fp);
        sscanf(line, "%s %s", NLGE, buf);
        printf("NGLE: %s, buf: %s\n", NLGE, buf);
        // fscanf(fp, "%c %s", NLGE, buf);
        if (!strcmp(NLGE, "N")) {
            printf("Objective read\n");
            strcpy(objective, buf);
        } else {
            constraint_signs[m] = NLGE;
            m++;
        }
        
        next_char = fgetc(fp);
        // next_char = fgetc(fp);
    }
    // printf("Next char is %c\n", next_char);
    command[0] = next_char;
    fgets(line, 100, fp);
    sscanf(line, "%s", &command[1]);
    if (strcmp(command, "COLUMNS")) {
        fprintf(stderr, "Wrong file format. Expected next command to be COLUMNS.\n");
        return 1;
    }
    
    next_char = fgetc(fp);
    printf("Next char is %c\n", next_char);
    c_float temp;
    
    char colchar[20], rowchar[20];

    while(next_char == ' ') {
        fgets(line, 100, fp);
        sscanf(line, "%s %s %le", colchar, rowchar, &temp);
        // fscanf(fp, "%c %s", NLGE, buf);
        if (!strcmp(rowchar, objective)) {
            printf("Objective read\n");
            
        } else {

        }
        
        next_char = fgetc(fp);
        // next_char = fgetc(fp);
    }

    // int character = fgetc(fp);
    // printf("The char was %c\n", character);
    
    fclose(fp);

    // // Load Q
    // fileno++;
    // fp = fopen(argv[fileno], "r");
    // if(fp == NULL) {
    //     fprintf(stderr, "Could not open file %s\n", argv[fileno]);
    // }
    // cholmod_sparse* Q = mtx_load_Q(fp, n);
    // fclose(fp);

    // // Load q
    // fileno++;
    // fp = fopen(argv[fileno], "r");
    // if(fp == NULL) {
    //     fprintf(stderr, "Could not open file %s\n", argv[fileno]);
    // }
    // c_float* q = mtx_load_dense(fp, n);
    // fclose(fp);

    // // Load lba
    // fileno++;
    // fp = fopen(argv[fileno], "r");
    // if(fp == NULL) {
    //     fprintf(stderr, "Could not open file %s\n", argv[fileno]);
    // }
    // c_float* bmin = mtx_load_dense(fp, m);
    // fclose(fp);

    // // Load uba
    // fileno++;
    // fp = fopen(argv[fileno], "r");
    // if(fp == NULL) {
    //     fprintf(stderr, "Could not open file %s\n", argv[fileno]);
    // }
    // c_float* bmax = mtx_load_dense(fp, m);
    // fclose(fp);



    // // Problem settings
    // QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));

    // // Structures
    // QPALMWorkspace *work; // Workspace
    // QPALMData *data;      // QPALMData

    // // Populate data
    // data    = (QPALMData *)c_malloc(sizeof(QPALMData));
    // data->n = n;
    // data->m = m;
    // data->q = q;
    // data->bmin = bmin;
    // data->bmax = bmax;
    // data->A = A;
    // data->Q = Q;

    // // Define Solver settings as default
    // qpalm_set_default_settings(settings);

    // // Setup workspace
    // cholmod_common c;
    // work = qpalm_setup(data, settings, &c);

    // // Solve Problem
    // qpalm_solve(work);

    // // printf("Solver status: %s\n", work->info->status);
    // // printf("Iter: %d\n", work->info->iter);
    // // printf("Iter_out: %d\n", work->info->iter_out);

    // // for (int i = 0; i < work->data->n; i++) {
    // //     printf("%f ", work->solution->x[i]);
    // // }
    // // printf("\n");

    // // #ifdef PROFILING
    // // printf("Setup time: %f\n", work->info->setup_time);
    // // printf("Solve time: %f\n", work->info->solve_time);
    // // printf("Run time: %f\n", work->info->run_time);
    // // #endif

    // // Clean workspace
    // CHOLMOD(start)(&work->chol->c);
    // CHOLMOD(free_sparse)(&data->Q, &c);
    // CHOLMOD(free_sparse)(&data->A, &c);
    // CHOLMOD(finish)(&work->chol->c);

    // qpalm_cleanup(work);
    // c_free(data->q);
    // c_free(data->bmin);
    // c_free(data->bmax);
    // c_free(data);
    // c_free(settings);
    return 0;
}