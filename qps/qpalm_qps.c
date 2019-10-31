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

    if (strcmp(command, check)) {
        fprintf(stderr, "Wrong command. Expected next command to be %s, but was %s.\n", check, command);
        return 1;
    }
    return 0;
}

long convert_to_long(char *s) {
    char c;
    long i, digit, number = 0;
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

    while(next_char == ' ') {
        fgets(line, 100, fp);
        sscanf(line, "%*s %s %le", rowchar, &temp);
        // row = convert_to_long(rowchar)-1;
        // switch (constraint_signs[row]) {
        //     case 'L':
        //         data->bmin[row] = data->bmax[row] - temp;
        //         break; 
        //     case 'G':
        //         data->bmax[row] = data->bmin[row] + temp;
        //         break;
        // }
        next_char = fgetc(fp);
    }

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
    Annz += n_bounds;

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

    // printf("Results: m = %lu, n = %lu, Qnnz = %lu, Annz = %lu\n", m, n, Qnnz, Annz);

    QPALMData* data = c_calloc(1, sizeof(QPALMData));
    data->m = m;
    data->n = n;
    data->c = 0;
    data->q = c_calloc(n, sizeof(c_float));
    size_t k;
    for (k = 0; k < n; k++) {
        data->q[k] = 0;
    }
    data->bmin = c_calloc(m, sizeof(c_float));
    for (k = 0; k < m; k++) {
        data->bmin[k] = 0;
    }
    data->bmax = c_calloc(m, sizeof(c_float));
    for (k = 0; k < m; k++) {
        data->bmax[k] = QPALM_INFTY;
    }

    cholmod_common c;
    CHOLMOD(start)(&c);
    data->A = CHOLMOD(allocate_sparse)(m, n, Annz, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
    data->Q = CHOLMOD(allocate_sparse)(n, n, Qnnz, TRUE, TRUE, -1, CHOLMOD_REAL, &c);

    c_float *Ax = data->A->x;
    c_int *Ai = data->A->i;
    c_int *Ap = data->A->p;

    c_float *Qx = data->Q->x;
    c_int *Qi = data->Q->i;
    c_int *Qp = data->Q->p;

    fp = fopen(argv[1], "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", argv[1]);
        return 1;
    }
    
    fgets(line, 100, fp);
    next_char = fgetc(fp);

    if(get_next_command_and_check(command, "ROWS", next_char, fp))
        return 1;

    next_char = fgetc(fp);

    char constraint_signs[m-n_bounds];
    size_t index = 0;
    while(next_char == ' ') {
        fgets(line, 100, fp);
        sscanf(line, "%s %s", NLGE, buf);
        if (strcmp(NLGE, "N")) {
            constraint_signs[index] = NLGE[0];
            index++;
        } 
        
        next_char = fgetc(fp);
    }

    if(get_next_command_and_check(command, "COLUMNS", next_char, fp))
        return 1;
    
    next_char = fgetc(fp);

    Ap[0] = 0;
    size_t elemA = 0;
    long row, col = 1, prev_col = 1;
    while(next_char == ' ') {
        fgets(line, 100, fp);
        sscanf(line, "%s %s %le %s %le", colchar, rowchar, &temp, second_rowchar, &temp2);
        row = convert_to_long(rowchar);
        col = convert_to_long(colchar);
        if (col > prev_col) {
            for (; prev_col < col; prev_col++) {
                if (prev_col <= n_bounds) { //take into account identity matrix for bounds
                    Ap[prev_col]++;
                    elemA++;
                }
                Ap[prev_col+1] = Ap[prev_col];
            }
        }

        if (!strcmp(rowchar, objective)) { 
            data->q[col-1] = temp;            
        } else {
            Ai[elemA] = row - 1;
            Ap[col]++;
            if (temp > QPALM_INFTY) {
                temp = QPALM_INFTY;
            } else if (temp < -QPALM_INFTY) {
                temp = -QPALM_INFTY;
            }
            Ax[elemA] = temp;
            elemA++;
        }

        if(!strcmp(second_rowchar, "")) {
        } else {
            row = convert_to_long(second_rowchar);
            Ai[elemA] = row - 1;
            Ap[col]++;
            if (temp2 > QPALM_INFTY) {
                temp2 = QPALM_INFTY;
            } else if (temp2 < -QPALM_INFTY) {
                temp2 = -QPALM_INFTY;
            }
            Ax[elemA] = temp2;
            elemA++;
            second_rowchar[0] = '\0';
        }
        
        next_char = fgetc(fp);
    }
    col = n;
    if (col > prev_col) {
        for (; prev_col < col; prev_col++) {
            if (prev_col <= n_bounds) { //take into account identity matrix for bounds
                Ap[prev_col]++;
                elemA++;
            }
            Ap[prev_col+1] = Ap[prev_col];
        }
    } else if ((col == prev_col) && prev_col <= n_bounds) { //take into account identity matrix for bounds
        Ap[prev_col]++;
    }

    
    if(get_next_command_and_check(command, "RHS", next_char, fp))
        return 1;

    next_char = fgetc(fp);
    k = 0;
    while(next_char == ' ') {
        fgets(line, 100, fp);
        sscanf(line, "%*s %s %le", rowchar, &temp);
        row = convert_to_long(rowchar)-1;
        if (!strcmp(rowchar, objective))
            data->c = -temp;
        else {
            switch (constraint_signs[row]) {
                case 'L':
                    data->bmax[row] = temp;
                    break; 
                case 'G':
                    data->bmin[row] = temp;
                    break;
                case 'E':
                    data->bmin[row] = temp;
                    data->bmax[row] = temp;
                    break;
            }
             
        }
        next_char = fgetc(fp);
    }

    long prev_row = row;
    if(get_next_command_and_check(command, "RANGES", next_char, fp))
        return 1;

    next_char = fgetc(fp);

    while(next_char == ' ') {
        fgets(line, 100, fp);
        sscanf(line, "%*s %s %le", rowchar, &temp);
        row = convert_to_long(rowchar)-1;
        switch (constraint_signs[row]) {
            case 'L':
                data->bmin[row] = data->bmax[row] - temp;
                break; 
            case 'G':
                data->bmax[row] = data->bmin[row] + temp;
                break;
        }
        next_char = fgetc(fp);
    }

    if(get_next_command_and_check(command, "BOUNDS", next_char, fp))
        return 1;

    next_char = fgetc(fp);
    k = prev_row; col = 1; long p; index = k; prev_row = -1;
    while(next_char == ' ') {
        fgets(line, 100, fp);
        sscanf(line, "%s %*s %s %le", bound_type, rowchar, &temp);
        row = convert_to_long(rowchar);
        if (!strcmp(bound_type, "UP")) {
            
            if(prev_row != row) {
                p = Ap[col];
                Ai[p-1] = row+k;
                Ax[p-1] = 1;
                col++;
                index++;
            }
            data->bmax[index] = temp;
        } else if (!strcmp(bound_type, "LO")) {
            
            if(prev_row != row) {
                p = Ap[col];
                Ai[p-1] = row+k;
                Ax[p-1] = 1;
                col++;
                index++;
            }
            data->bmin[index] = temp;
        }
        prev_row = row;
        next_char = fgetc(fp);
    }

    if(get_next_command_and_check(command, "QUADOBJ", next_char, fp))
        return 1;
    next_char = fgetc(fp);
    prev_col = 1;
    size_t elemQ = 0;
    Qp[0] = 0;
    while(next_char == ' ') {
        fgets(line, 100, fp);
        sscanf(line, "%s %s %le", colchar, rowchar, &temp);
        col = convert_to_long(colchar);
        row = convert_to_long(rowchar);

        Qi[elemQ] = --row;
        
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
        Qx[elemQ] = temp;

        elemQ++;
        next_char = fgetc(fp);
    }


    if(get_next_command_and_check(command, "ENDATA", next_char, fp))
        return 1;

    fclose(fp);

    CHOLMOD(finish)(&c);
    

    // Setup and solve problem

    // Problem settings
    QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    // Define Solver settings as default
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->max_iter = 10000;
    // settings->verbose = FALSE;

    QPALMWorkspace *work = qpalm_setup(data, settings, &c);

    // Solve Problem
    qpalm_solve(work);

    strcpy(name, &(argv[1][78]));
    
    // name = &name[78];
    name[strlen(name)-4] = '\0';
    // printf("%s\n", name);
    fp = fopen("out.tex", "a");
    fprintf(fp, "%s & %lu & %lu & %lu & %lu & %lu & %le \\\\", name, n, m, Annz, Qnnz, work->info->iter, work->info->objective);
    fclose(fp);

    // Clean workspace
    CHOLMOD(start)(&c);
    CHOLMOD(free_sparse)(&data->Q, &c);
    CHOLMOD(free_sparse)(&data->A, &c);
    CHOLMOD(finish)(&c);
    qpalm_cleanup(work);

    c_free(data->q);
    c_free(data->bmin);
    c_free(data->bmax);
    c_free(data);
    c_free(settings);

    return 0;
}