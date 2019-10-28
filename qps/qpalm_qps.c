#include "qpalm.h"
#include "constants.h"
#include "global_opts.h"
#include "cholmod.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int get_next_command_and_check(char* command, char* check, char next_char, FILE* fp) {
    command[0] = next_char;
    char line[100];
    fgets(line, 100, fp);
    sscanf(line, "%s", &command[1]);
    // printf("Next command: %s\n", command);

    if (strcmp(command, check)) {
        fprintf(stderr, "Wrong command. Expected next command to be %s, but was %s.\n", check, command);
        return 1;
    }
    return 0;
}

long convert_to_long(char *s) {
    char c;
    long i, digit, number;
    for(i=0;i<strlen(s);i++)
    {
    	c = s[i];
    	if(c>='0' && c<='9') //to confirm it's a digit
    	{
    		digit = c - '0';
    		number = number*10 + digit;
    	}
    }
    return number;
}

int main(int argc, char*argv[]){

    if (argc != 2) {
        fprintf(stderr, "Wrong number of arguments. Correct usage is qpalm_qps problem.qps\n");
    }
    
    // Load problem data
    size_t n, m, Annz, Qnnz, n_bounds;
    n = 0; m = 0; Annz = 0; Qnnz = 0; n_bounds = 0;
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
    }

    if (strcmp(command, "NAME")) {
        fprintf(stderr, "Wrong file format. Expected first line to contain NAME problem_name nnz.\n");
        return 1;
    }

    printf("Reading problem %s\n", name);

    char next_char;
    next_char = fgetc(fp);

    if(get_next_command_and_check(command, "ROWS", next_char, fp))
        return 1;

    next_char = fgetc(fp);
    char NLGE[1], buf[20], objective[20];

    while(next_char == ' ') {
        fgets(line, 100, fp);
        sscanf(line, "%s %s", NLGE, buf);
        if (!strcmp(NLGE, "N")) {
            strcpy(objective, buf);
        } else {
            m++;
        }
        
        next_char = fgetc(fp);
    }
    if(get_next_command_and_check(command, "COLUMNS", next_char, fp))
        return 1;
    
    next_char = fgetc(fp);
    
    c_float temp, temp2;
    char colchar[20], rowchar[20], prev_colchar[20], second_rowchar[20];
    prev_colchar[0] = '\0';
    second_rowchar[0] = '\0';

    while(next_char == ' ') {
        fgets(line, 100, fp);
        sscanf(line, "%s %s %le %s %le", colchar, rowchar, &temp, second_rowchar, &temp2);

        if (strcmp(colchar, prev_colchar)) {
            n++;
            strcpy(prev_colchar, colchar);
        }

        if (!strcmp(rowchar, objective)) {            
        } else {
            Annz++;
        }

        if (!strcmp(second_rowchar, objective)){
        } else if(!strcmp(second_rowchar, "")) {
        } else {
            Annz++;
            second_rowchar[0] = '\0';
        }
        
        next_char = fgetc(fp);
    }
    if(get_next_command_and_check(command, "RHS", next_char, fp))
        return 1;

    next_char = fgetc(fp);
    while(next_char == ' ') {
        fgets(line, 100, fp);
        next_char = fgetc(fp);
    }

    if(get_next_command_and_check(command, "RANGES", next_char, fp))
        return 1;

    next_char = fgetc(fp);
    if(get_next_command_and_check(command, "BOUNDS", next_char, fp))
        return 1;
    next_char = fgetc(fp);

    char bound_type[20];
    char prev_rowchar[20];
    prev_rowchar[0] = '\0';
    while(next_char == ' ') {
        fgets(line, 100, fp);
        sscanf(line, "%s %*s %s %le", bound_type, rowchar, &temp);
        if (strcmp(bound_type, "FR") && strcmp(rowchar, prev_rowchar)) {
            m++;
            n_bounds++;
            strcpy(prev_rowchar, rowchar);
        }
        next_char = fgetc(fp);
    }

    if(get_next_command_and_check(command, "QUADOBJ", next_char, fp))
        return 1;
    next_char = fgetc(fp);
    while(next_char == ' ') {
        fgets(line, 100, fp);
        next_char = fgetc(fp);
        Qnnz++;
    }


    if(get_next_command_and_check(command, "ENDATA", next_char, fp))
        return 1;
    fclose(fp);

    printf("Results: m = %lu, n = %lu, Qnnz = %lu, Annz = %lu\n", m, n, Qnnz, Annz);

    c_float *q = c_calloc(n, sizeof(c_float));
    c_float *bmin = c_calloc(m, sizeof(c_float));
    if (n_bounds > 0) {
        for (size_t k = m-n_bounds; k < m; k++) {
            bmin[k] = 0;
        }
    }
    c_float *bmax = c_calloc(m, sizeof(c_float));
    

    fp = fopen(argv[1], "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", argv[1]);
        return 1;
    }
    
    cholmod_common c;
    CHOLMOD(start)(&c);

    cholmod_sparse *A = CHOLMOD(allocate_sparse)(m, n, Annz, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    cholmod_sparse *Q = CHOLMOD(allocate_sparse)(n, n, Qnnz, TRUE, TRUE, -1, CHOLMOD_REAL, &c);



    CHOLMOD(finish)(&c);
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