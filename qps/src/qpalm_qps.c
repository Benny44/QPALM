#include "qpalm.h"
#include "constants.h"
#include "global_opts.h"
#include "cholmod.h"
#include "qpalm_qps.h"
#include "qps_conversion.h"
#include "index_hash.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Print a cholmod matrix so the output can be entered into matlab */
void print_cholmod_matlab(cholmod_sparse *M) {
    printf("M = sparse(%ld, %ld);", M->nrow, M->ncol);
    size_t col, index = 0, row;
    double *Mx = M->x;
    long int *Mi = M->i;
    long int *Mp = M->p;


    for (col = 1; col <= M->ncol; col++) {
        for (row = Mp[col-1]; row < Mp[col]; row++) {
            printf("M(%ld, %ld) = %.16le;", Mi[index]+1, col, Mx[index]);
            index++;
        }
    }
    printf("\n");
}

void print_dense_vector_matlab(double* x, size_t len) {
    size_t k;
    printf("x = zeros(%lu, 1);", len);
    for (k = 0; k < len; k++) {
        printf("x(%lu) = %.16le;", k+1, x[k]);
    }
    printf("\n");
}

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

int get_next_command(char* command, char next_char, FILE* fp) {
    command[0] = next_char;
    char line[100];
    fgets(line, 100, fp);
    sscanf(line, "%s", &command[1]);

    return strcmp(command, "ENDATA");
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

int get_sizes_and_check_format(FILE *fp, QPALMData *data, struct index_table **free_bounds) {

    size_t n, m, Annz, Qnnz, n_bounds;
    n = 0; m = 0; Annz = 0; Qnnz = 0; n_bounds = 0;

    char next_char;
    char line[100], command[20];
    int old_format_detected = FALSE;
    // char *file_copy = NULL;
    char NLGE[1], buf[20], buf2[20], objective[20];
    buf2[0] = '\0';
    c_float temp, temp2;
    char colchar[20], rowchar[20], prev_colchar[20], second_rowchar[20];
    prev_colchar[0] = '\0';
    second_rowchar[0] = '\0';

    // c_int *bounds = c_calloc(n, sizeof(c_int));
    size_t k;
    
    char bound_type[20];
    char prev_rowchar[20];
    prev_rowchar[0] = '\0';
    long row;
    

    next_char = fgetc(fp);
    /*First pass through the file to get the sizes*/
    while(get_next_command(command, next_char, fp)){
        next_char = fgetc(fp);
        if (!strcmp(command, "ROWS")) {
    
            while(next_char == ' ') {
                fgets(line, 100, fp);
                sscanf(line, "%s %s %s", NLGE, buf, buf2);
                if (strcmp(buf2, "")){
                    printf("Old qps format detected. First performing conversion to new format.\n");
                    old_format_detected = TRUE;
                    break;
                }
                if (!strcmp(NLGE, "N")) {
                    strcpy(objective, buf);
                } else {
                    m++;
                }
                
                next_char = fgetc(fp);
            }
            if (old_format_detected) {
                break;
            }

        } else if (!strcmp(command, "COLUMNS")) {

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
            n_bounds = n;
            Annz += n_bounds;
            m += n_bounds;
            *free_bounds = create_index_table(c_max(n/5,1));

        } else if (!strcmp(command, "RHS") || !strcmp(command, "RANGES")) {
            while(next_char == ' ') {
                fgets(line, 100, fp);
                next_char = fgetc(fp);
            }

        } else if (!strcmp(command, "BOUNDS")) {
            
            while(next_char == ' ') {
                fgets(line, 100, fp);
                sscanf(line, "%s %*s %s %le", bound_type, colchar, &temp);
                if (!strcmp(bound_type, "FR")) {
                    // col = convert_to_long(colchar)-1;
                    // bounds[col] = FALSE;
                    insert(*free_bounds, colchar, 0, ' ');
                    n_bounds--;
                    m--;
                    Annz--;
                } 

                next_char = fgetc(fp);
            }
            

        } else if (!strcmp(command, "QUADOBJ")) {
            while(next_char == ' ') {
                fgets(line, 100, fp);
                next_char = fgetc(fp);
                Qnnz++;
            }
        }
    }

    fclose(fp);
    if (!old_format_detected) {
        data->m = m;
        data->n = n;
        data->c = 0;
        data->q = c_calloc(n, sizeof(c_float));
        
        for (k = 0; k < n; k++) {
            data->q[k] = 0;
        }
        data->bmin = c_calloc(m, sizeof(c_float));
        data->bmax = c_calloc(m, sizeof(c_float));    

        cholmod_common c;
        CHOLMOD(start)(&c);
        data->A = CHOLMOD(allocate_sparse)(m, n, Annz, TRUE, TRUE, 0, CHOLMOD_REAL, &c);
        data->Q = CHOLMOD(allocate_sparse)(n, n, Qnnz, TRUE, TRUE, -1, CHOLMOD_REAL, &c);
        CHOLMOD(finish)(&c);
    }
    return old_format_detected;
}

void read_data(FILE* fp, QPALMData *data, struct index_table* row_index_table,
                struct index_table* col_index_table, struct index_table* free_bounds) {

    char next_char;
    char NLGE[1], buf[20], buf2[20], objective[20];
    buf2[0] = '\0';
    objective[0] = '\0';
    c_float temp, temp2;
    char colchar[20], rowchar[20], prev_colchar[20], second_rowchar[20];
    prev_colchar[0] = '\0';
    second_rowchar[0] = '\0';

    size_t k;
    
    char bound_type[20];
    char prev_rowchar[20];
    prev_rowchar[0] = '\0';
    c_int row;

    char line[100], command[20], name[50];
    size_t n, m, n_bounds;
    m = data->m;
    n = data->n;
    n_bounds = n - length_table(free_bounds, c_max(n/5,1));

    c_float *Ax = data->A->x;
    c_int *Ai = data->A->i;
    c_int *Ap = data->A->p;

    c_float *Qx = data->Q->x;
    c_int *Qi = data->Q->i;
    c_int *Qp = data->Q->p;

    size_t index = 0;
    long bounds_row = m-n_bounds, col, prev_col = 1;

    char constraint_sign;
    struct node *row_node;
    struct node *col_node;

    fgets(line, 100, fp);
    next_char = fgetc(fp);
    while (get_next_command(command, next_char, fp)) {
        next_char = fgetc(fp);
        if (!strcmp(command, "ROWS")) {
            index = 0;
            while(next_char == ' ') {
                fgets(line, 100, fp);
                sscanf(line, "%s %s", NLGE, buf);
                constraint_sign = NLGE[0];
                if (constraint_sign == 'N') {
                    strcpy(objective, buf);
                } else {
                    insert(row_index_table, buf, index, constraint_sign);
                    switch (constraint_sign) {
                        case 'L':
                            data->bmax[index] = 0;
                            data->bmin[index] = -QPALM_INFTY;
                            break; 
                        case 'G':
                            data->bmin[index] = 0;
                            data->bmax[index] = QPALM_INFTY;
                            break;
                        case 'E':
                            data->bmin[index] = 0;
                            data->bmax[index] = 0;
                            break;
                    
                    }
                    index++;
                } 
                next_char = fgetc(fp);
            }

            for (k = m-n_bounds; k < m; k++) {
                    data->bmin[k] = 0;
                    data->bmax[k] = QPALM_INFTY;
            }
        } else if (!strcmp(command, "COLUMNS")) {
            Ap[0] = 0;
            size_t elemA = 0;
            prev_col = 1;
            col = 0;
            prev_colchar[0] = '\0';

            while(next_char == ' ') {
                fgets(line, 100, fp);
                sscanf(line, "%s %s %le %s %le", colchar, rowchar, &temp, second_rowchar, &temp2);

                if (strcmp(colchar, prev_colchar)) {
                    col++;
                    insert(col_index_table, colchar, col, ' ');
                    for (; prev_col < col; prev_col++) {
                        if (lookup(free_bounds, prev_colchar) == NULL) {
                            Ai[Ap[prev_col]] = bounds_row;
                            Ax[Ap[prev_col]] = 1;
                            Ap[prev_col]++;
                            elemA++;
                            bounds_row++;
                        }
                        Ap[prev_col+1] = Ap[prev_col];
                    }
                    strcpy(prev_colchar, colchar);
                }

                if (!strcmp(rowchar, objective)) { 
                    data->q[col-1] = temp;            
                } else {
                    row_node = lookup(row_index_table, rowchar);
                    if (row_node == NULL) {
                        printf("Line: %s\n", line);
                        printf("Rowchar: %s\n", rowchar);
                        printf("Objective: %s\n", objective);
                    }
                    row = row_node->index;
                    Ai[elemA] = row;
                    Ap[col]++;
                    if (temp > QPALM_INFTY) {
                        temp = QPALM_INFTY;
                    } else if (temp < -QPALM_INFTY) {
                        temp = -QPALM_INFTY;
                    }
                    Ax[elemA] = temp;
                    elemA++;
                }

                if(strcmp(second_rowchar, "") && strcmp(second_rowchar, objective)) {
                    // printf("Second row: %s\n", second_rowchar);
                    row_node = lookup(row_index_table, second_rowchar);
                    // if (row_node == NULL) printf("Row not found\n");
                    row = row_node->index;
                    // row = convert_to_long(second_rowchar)-1;
                    Ai[elemA] = row;
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
            if ((col == prev_col) && lookup(free_bounds, prev_colchar) == NULL) { //take into account identity matrix for bounds
                Ai[Ap[prev_col]] = bounds_row;
                Ax[Ap[prev_col]] = 1;
                Ap[prev_col]++;
                elemA++;
                bounds_row++;
            }
        } else if (!strcmp(command, "RHS")) {
            // printf("RHS\n");
            while(next_char == ' ') {
                fgets(line, 100, fp);
                sscanf(line, "%*s %s %le %s %le", rowchar, &temp, second_rowchar, &temp2);
                // row = convert_to_long(rowchar)-1;

                if (!strcmp(rowchar, objective))
                    data->c = -temp;
                else {
                    row_node = lookup(row_index_table, rowchar);
                    row = row_node->index;
                    
                    switch (row_node->sign) {
                        case 'L':
                            data->bmax[row] = temp;
                            data->bmin[row] = -QPALM_INFTY;
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

                if (!strcmp(second_rowchar, objective))
                    data->c = -temp2;
                else if (strcmp(second_rowchar, "")){
                    row_node = lookup(row_index_table, second_rowchar);
                    row = row_node->index;
                    
                    switch (row_node->sign) {
                        case 'L':
                            data->bmax[row] = temp2;
                            data->bmin[row] = -QPALM_INFTY;
                            break; 
                        case 'G':
                            data->bmin[row] = temp2;
                            break;
                        case 'E':
                            data->bmin[row] = temp2;
                            data->bmax[row] = temp2;
                            break;
                    }
                }

                next_char = fgetc(fp);
            }
        } else if (!strcmp(command, "RANGES")) {
                        // printf("RANGES\n");

            while(next_char == ' ') {
                fgets(line, 100, fp);
                sscanf(line, "%*s %s %le %s %le", rowchar, &temp, second_rowchar, &temp2);
                // row = convert_to_long(rowchar)-1;
                row_node = lookup(row_index_table, rowchar);
                row = row_node->index;
                switch (row_node->sign) {
                    case 'L':
                        data->bmin[row] = data->bmax[row] - temp;
                        break; 
                    case 'G':
                        data->bmax[row] = data->bmin[row] + temp;
                        break;
                }

                if (strcmp(second_rowchar, "")) {
                    row_node = lookup(row_index_table, second_rowchar);
                    row = row_node->index;
                    switch (row_node->sign) {
                        case 'L':
                            data->bmin[row] = data->bmax[row] - temp2;
                            break; 
                        case 'G':
                            data->bmax[row] = data->bmin[row] + temp2;
                            break;
                    }
                }

                next_char = fgetc(fp);
            }
        } else if (!strcmp(command, "BOUNDS")) {

            index = m-n_bounds;
            // print_table(col_index_table, c_max(n/5, 1));
            while(next_char == ' ') {
                fgets(line, 100, fp);
                sscanf(line, "%s %*s %s %le", bound_type, colchar, &temp);
                // col = convert_to_long(colchar)-1;
                col_node = lookup(col_index_table, colchar);
                col = col_node->index - 1;
                
                // printf("Found %s with index %ld\n", colchar, col_node->index);

                if (!strcmp(bound_type, "FR")) {
                    index--;
                } else if (!strcmp(bound_type, "UP")) {
                    data->bmax[index+col] = temp;
                } else if (!strcmp(bound_type, "LO")) {
                    // printf("Inserting in bmin at %ld + %ld\n", index, col);
                    data->bmin[index+col] = temp;
                } else if (!strcmp(bound_type, "FX")) {
                    data->bmin[index+col] = temp;
                    data->bmax[index+col] = temp;
                }
                
                next_char = fgetc(fp);
            }
        } else if (!strcmp(command, "QUADOBJ")) {
            prev_col = 1;
            size_t elemQ = 0;
            Qp[0] = 0;
            while(next_char == ' ') {
                fgets(line, 100, fp);
                sscanf(line, "%s %s %le", colchar, rowchar, &temp);
                // col = convert_to_long(colchar);
                // row = convert_to_long(rowchar);
                col_node = lookup(col_index_table, colchar);
                col = col_node->index;
                row_node = lookup(col_index_table, rowchar);
                row = row_node->index;

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
            col = n;
            if (col > prev_col) {
                    for (; prev_col < col; prev_col++) {
                        Qp[prev_col+1] = Qp[prev_col];
                    }          
            }
        }
    }
}

void print_out_latex(QPALMData *data, QPALMInfo *info, char* file){
    char line[100], name[20];
    strcpy(line, file);
    size_t last_slash = 0;
    size_t k;
    for (k = 0; k < strlen(line); k++) {
        if (line[k] == '/') /*assume linux directory*/
            last_slash = k+1;
    }
    strcpy(name, &line[last_slash]);
    name[strlen(name)-4] = '\0'; /*delete .qps*/
    FILE *fp = fopen("out.tex", "a");
    fprintf(fp, "%s & %lu & %lu & %lu & %lu & %lu & %le \\\\ \n", name, data->n, data->m, data->A->nzmax, data->Q->nzmax, info->iter, info->objective);
    fclose(fp);
}

void convert_underscore_to_dash(char* name) {
    char c;
    while ((c = *name) != '\0') {
        if (c == '_')
            *name = '-';
        name++;
    }
}

void print_out_bpmpd(QPALMData *data, QPALMInfo *info, char* file){
    char line[100], name[20], buf[20];
    c_float temp = -QPALM_INFTY;
    strcpy(line, file);
    size_t last_slash = 0;
    size_t k;
    for (k = 0; k < strlen(line); k++) {
        if (line[k] == '/') /*assume linux directory*/
            last_slash = k+1;
    }
    strcpy(name, &line[last_slash]);
    name[strlen(name)-4] = '\0'; /*delete .qps*/
    FILE *fp = fopen("out.tex", "a");

    /*Get the objective from BPMPD*/
    FILE *fp_bpmpd = fopen("/home/ben/Documents/Projects/QPALM/qps/BPMPD.txt", "r");
    if (fp_bpmpd == NULL) {
        printf("Could not open\n");
        return;
    }

    fgets(line, 100, fp_bpmpd);
    sscanf(line, "%s %le", buf, &temp);
    while (strcmp(buf, name)) {
        if (fgets(line, 100, fp_bpmpd) != NULL) {
            sscanf(line, "%s %le", buf, &temp);
            // printf("Buf: %s, temp: %le\n", buf, temp);
        }
        else
            break;
    }

    if (!strcmp(buf, name)) {
        convert_underscore_to_dash(name);
        fprintf(fp, "%s & %lu & %lu & %lu & %lu & %lu & %le & %le \\\\ \n", name, data->n, data->m, data->A->nzmax, data->Q->nzmax, info->iter, info->objective, temp);
        // fprintf(fp, "%s & %le \\\\ \n", name, temp);
    } else
        printf("Name not found in BPMPD file: %s\n");

    // fprintf(fp, "%s & %lu & %lu & %lu & %lu & %lu & %le \\\\ \n", name, data->n, data->m, data->A->nzmax, data->Q->nzmax, info->iter, info->objective);
    fclose(fp);
    fclose(fp_bpmpd);
}

int main(int argc, char*argv[]){

    if (argc != 2) {
        fprintf(stderr, "Wrong number of arguments. Correct usage is qpalm_qps problem.qps\n");
    }


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
    char *file_copy = NULL;    
    struct index_table * free_bounds;
    QPALMData* data = c_calloc(1, sizeof(QPALMData));

    /* First pass through the file to get the sizes. If an old QPS-format is 
       detected (which allows for spaces in names), first a conversion is 
       performed to the new format and then the (new) file is read again.
    */
    if (get_sizes_and_check_format(fp, data, &free_bounds)) {
        file_copy = convert_qps_to_new_format(argv[1]);
        fp = fopen(file_copy, "r");
        if(fp == NULL) {
            fprintf(stderr, "Could not open file %s\n", file_copy);
            return 1;
        }
        fgets(line, 100, fp);
        get_sizes_and_check_format(fp, data, &free_bounds);
    }

    
    /* Go through file a second time to read in the data*/
    if (file_copy == NULL) {
        fp = fopen(argv[1], "r");
        if(fp == NULL) {
            fprintf(stderr, "Could not open file %s\n", argv[1]);
            return 1;
        }
    } else {
        fp = fopen(file_copy, "r");
        if(fp == NULL) {
            fprintf(stderr, "Could not open file %s\n", file_copy);
            return 1;
        }
    } 

    size_t m, n, n_bounds;
    m = data->m;
    n = data->n;
    n_bounds = n - length_table(free_bounds, c_max(n/5,1));

    struct index_table* row_index_table = create_index_table(c_max((m-n_bounds)/5, 1));
    struct index_table* col_index_table = create_index_table(c_max(n/5, 1));
    read_data(fp, data, row_index_table, col_index_table, free_bounds);

    fclose(fp);
    if (file_copy) c_free(file_copy);

    printf("Reading successful.\n");

    // print_cholmod_matlab(data->Q);
    // print_cholmod_matlab(data->A);
    // print_dense_vector_matlab(data->bmin, m);
    // print_dense_vector_matlab(data->bmax, m);


    /* Setup problem */
    QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->eps_dual_inf = 1e-6;
    settings->eps_prim_inf = 1e-6;
    settings->max_iter = 10000;
    settings->verbose = TRUE;
    // settings->scaling = 2;
    // settings->proximal = TRUE;

    cholmod_common c;
    QPALMWorkspace *work = qpalm_setup(data, settings, &c);

    /* Solve Problem */
    qpalm_solve(work);
    
    // printf("Iter: %ld\n", work->info->iter);
    // printf("Status: %s\n", work->info->status);
    // printf("Objective: %le\n", work->info->objective);

    // print_out_latex(data, work->info, argv[1]);
    print_out_bpmpd(data, work->info, argv[1]);

    // // Clean workspace
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

    // c_free(bounds);
    free_index_table(row_index_table, c_max((m-n_bounds)/5,1));
    free_index_table(col_index_table, c_max((n)/5,1));
    free_index_table(free_bounds, c_max((n)/5,1));

    return 0;
}