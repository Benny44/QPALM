#include "qpalm.h"
#include "constants.h"
#include "global_opts.h"
#include "cholmod.h"
#include "qpalm_qps.h"
#include "qps_conversion.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


struct node{
    char* key; /*name of the row/col*/
    c_int index; /*index of the row/col */
    char sign; /*constraint sign for a row*/
    struct node *next;
};
struct index_table{
    c_int size;
    struct node **list;
};
struct index_table *create_index_table(c_int size){
    struct index_table *t = (struct table*)malloc(sizeof(struct index_table));
    t->size = size;
    t->list = (struct node**)malloc(sizeof(struct node*)*size);
    c_int i;
    for(i=0; i<size; i++)
        t->list[i] = NULL;
    return t;
}
c_int hashcode(struct index_table *t,char* key){
    c_int code = 0;
    while(*key) code += (c_int)*key++;
    return code % t->size;
}
void insert(struct index_table *t, char* key, c_int index, char sign){
    c_int pos = hashcode(t,key);
    struct node *list = t->list[pos];
    struct node *new = (struct node*)malloc(sizeof(struct node));
    struct node *temp = list;
    while (temp) temp = temp->next;
    new->key = (char*)malloc(sizeof(char)*(strlen(key)+1));
    strcpy(new->key, key);
    new->index = index;
    new->sign = sign;
    new->next = list;
    t->list[pos] = new;
}
struct node* lookup(struct index_table *t,char* key){
    c_int pos = hashcode(t,key);
    struct node *list = t->list[pos];
    struct node *temp = list;
    while(temp){
        if (strcmp(key, temp->key)) temp=temp->next;
        else return temp;
    }
    return NULL;
}
void print_table(struct index_table *t, c_int size){
    struct node *temp;
    struct node *list;
    c_int pos;
    for (pos = 0; pos < size; pos++) {
        list = t->list[pos];
        temp = list;
        while (temp) {
            printf("Entry: %s, %ld\n", temp->key, temp->index);
            temp = temp->next;
        }
    }
}

void free_index_table(struct index_table *t, c_int size) {
    struct node *temp, *temp_next;
    struct node *list;
    c_int pos;
    for (pos = 0; pos < size; pos++) {
        list = t->list[pos];
        temp = list;
        while (temp) {
            temp_next = temp->next;
            c_free(temp->key);
            c_free(temp);
            temp = temp_next;
            // printf("Entry: %s, %ld\n", temp->key, temp->index);
            // temp = temp->next;
        }
    }
    c_free(t->list);
    c_free(t);
}
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
    
    
    int old_format_detected = FALSE;
    char *file_copy = NULL;
    char NLGE[1], buf[20], buf2[20], objective[20];
    buf2[0] = '\0';
    c_float temp, temp2;
    char colchar[20], rowchar[20], prev_colchar[20], second_rowchar[20];
    prev_colchar[0] = '\0';
    second_rowchar[0] = '\0';

    // c_int *bounds = c_calloc(n, sizeof(c_int));
    c_int *bounds;
    size_t k;
    
    char bound_type[20];
    char prev_rowchar[20];
    prev_rowchar[0] = '\0';
    long row;
    
    struct index_table * free_bounds;

    /*First pass through the file to get the sizes*/
    while(get_next_command(command, next_char, fp)){
        next_char = fgetc(fp);

        if (!strcmp(command, "ROWS")) {
    
            while(next_char == ' ') {
                fgets(line, 100, fp);
                sscanf(line, "%s %s %s", NLGE, buf, buf2);
                if (strcmp(buf2, "")){
                    // printf("Old qps format detected. First performing conversion to new format.\n");
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
            free_bounds = create_index_table(c_max(n/5,1));

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
                    insert(free_bounds, colchar, 0, ' ');
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
    
    if (old_format_detected) {
        m = 0;
        buf2[0] = '\0';
        file_copy = convert_qps_to_new_format(argv[1]);
        fp = fopen(file_copy, "r");
        if(fp == NULL) {
            fprintf(stderr, "Could not open file %s\n", file_copy);
            return 1;
        }
        fgets(line, 100, fp);
        next_char = fgetc(fp);
        old_format_detected = FALSE;

        /*First pass through the file to get the sizes*/
        while(get_next_command(command, next_char, fp)){
            next_char = fgetc(fp);

            if (!strcmp(command, "ROWS")) {
        
                while(next_char == ' ') {
                    fgets(line, 100, fp);
                    sscanf(line, "%s %s %s", NLGE, buf, buf2);
                    if (strcmp(buf2, "")){
                        // printf("Old qps format detected. First performing conversion to new format.\n");
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
                free_bounds = create_index_table(c_max(n/5,1));

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
                        insert(free_bounds, colchar, 0, ' ');
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
    }

    // printf("Results: m = %lu, n = %lu, Qnnz = %lu, Annz = %lu\n", m, n, Qnnz, Annz);

    struct index_table* row_index_table = create_index_table(c_max((m-n_bounds)/5, 1));
    struct index_table* col_index_table = create_index_table(c_max(n/5, 1));

    QPALMData* data = c_calloc(1, sizeof(QPALMData));
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

    c_float *Ax = data->A->x;
    c_int *Ai = data->A->i;
    c_int *Ap = data->A->p;

    c_float *Qx = data->Q->x;
    c_int *Qi = data->Q->i;
    c_int *Qp = data->Q->p;

    if (file_copy == NULL) fp = fopen(argv[1], "r");
    else fp = fopen(file_copy, "r");
    if(fp == NULL) {
        fprintf(stderr, "Could not open file %s\n", argv[1]);
        return 1;
    }

    // char constraint_signs[m-n_bounds];
    size_t index = 0;
    long bounds_row = m-n_bounds, col, prev_col = 1;
    
    fgets(line, 100, fp);
    next_char = fgetc(fp);

    char constraint_sign;
    struct node *row_node;
    struct node *col_node;

    // print_table(free_bounds, c_max(n/5,1));

    while (get_next_command(command, next_char, fp)) {
        next_char = fgetc(fp);
        if (!strcmp(command, "ROWS")) {
            index = 0;
            while(next_char == ' ') {
                fgets(line, 100, fp);
                sscanf(line, "%s %s", NLGE, buf);
                constraint_sign = NLGE[0];
                if (constraint_sign == 'E') { //default in equality for bmax = bmin = 0
                    data->bmax[index] = 0;
                }
                if (constraint_sign != 'N') {
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
                    // constraint_signs[index] = NLGE[0];
                    index++;
                } 
                next_char = fgetc(fp);
            }

            for (k = m-n_bounds; k < m; k++) {
                // if (k >= m-n_bounds) {
                    data->bmin[k] = 0;
                    data->bmax[k] = QPALM_INFTY;
                // }
                // else switch (constraint_signs[k]) {
                //         case 'L':
                //             data->bmax[k] = 0;
                //             data->bmin[k] = -QPALM_INFTY;
                //             break; 
                //         case 'G':
                //             data->bmin[k] = 0;
                //             data->bmax[k] = QPALM_INFTY;
                //             break;
                //         case 'E':
                //             data->bmin[k] = 0;
                //             data->bmax[k] = 0;
                //             break;
                    
                // }
            }
        } else if (!strcmp(command, "COLUMNS")) {
            // print_table(row_index_table, c_max((m-n_bounds)/5,1));
            Ap[0] = 0;
            size_t elemA = 0;
            prev_col = 1;
            col = 0;
            prev_colchar[0] = '\0';

            while(next_char == ' ') {
                fgets(line, 100, fp);
                sscanf(line, "%s %s %le %s %le", colchar, rowchar, &temp, second_rowchar, &temp2);
                // row = convert_to_long(rowchar)-1;
                
                // col_node = lookup(col_index_table, colchar);
                // col = col_node->index;
                // col = convert_to_long(colchar);

                if (strcmp(colchar, prev_colchar)) {
                    // n++;
                    col++;
                    // printf("Inserting: %s as column %ld \n", colchar, col);
                    insert(col_index_table, colchar, col, ' ');
                    for (; prev_col < col; prev_col++) {
                        // if (bounds[prev_col-1]) { //take into account identity matrix for bounds
                        if (lookup(free_bounds, prev_colchar) == NULL) {
                        // printf("Prev_col: %ld, bounds = %ld\n", prev_col, bounds[prev_col-1]);
                        // p = Ap[col+1];
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

                // if (col > prev_col) {
                    
                // }

                if (!strcmp(rowchar, objective)) { 
                    data->q[col-1] = temp;            
                } else {
                    row_node = lookup(row_index_table, rowchar);
                    if (row_node == NULL) {
                        printf("Line: %s\n", line);
                        printf("Rowchar: %s\n", rowchar);
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
            // col = n;
            // if (col > prev_col) {
            //     for (; prev_col < col; prev_col++) {
            //         if (bounds[prev_col-1]) { //take into account identity matrix for bounds
            //             Ai[Ap[prev_col]] = bounds_row;
            //                 Ax[Ap[prev_col]] = 1;
            //                 Ap[prev_col]++;
            //                 elemA++;
            //                 bounds_row++;
            //         }
            //         Ap[prev_col+1] = Ap[prev_col];
            //     }
            // } else 
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

    fclose(fp);
    if (file_copy) c_free(file_copy);

    CHOLMOD(finish)(&c);

    // print_cholmod_matlab(data->Q);
    // print_cholmod_matlab(data->A);
    // print_dense_vector_matlab(data->bmin, m);
    // print_dense_vector_matlab(data->bmax, m);


    // Problem settings
    QPALMSettings *settings = (QPALMSettings *)c_malloc(sizeof(QPALMSettings));
    // Define Solver settings as default
    qpalm_set_default_settings(settings);
    settings->eps_abs = 1e-6;
    settings->eps_rel = 1e-6;
    settings->eps_dual_inf = 1e-6;
    settings->eps_prim_inf = 1e-6;
    settings->max_iter = 1000;
    settings->verbose = FALSE;
    settings->scaling = FALSE;
    // settings->proximal = TRUE;

    QPALMWorkspace *work = qpalm_setup(data, settings, &c);

    // printf("Before solve\n");

    // Solve Problem
    qpalm_solve(work);
    // printf("After solve\n");

    // for (k = 0; k < n; k++) {
    //     printf("x_sol[%ld] = %le\n", k, work->solution->x[k]);
    // }
    printf("Iter: %ld\n", work->info->iter);
    printf("Status: %s\n", work->info->status);
    printf("Objective: %le\n", work->info->objective);
    
    strcpy(line, argv[1]);
    size_t last_slash = 0;
    for (k = 0; k < strlen(line); k++) {
        if (line[k] == '/') /*assume linux directory*/
            last_slash = k+1;
    }
    strcpy(name, &line[last_slash]);
    // name = &name[78];
    name[strlen(name)-4] = '\0'; /*delete .qps*/
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

    // c_free(bounds);
    free_index_table(row_index_table, c_max((m-n_bounds)/5,1));
    free_index_table(col_index_table, c_max((n)/5,1));
    free_index_table(free_bounds, c_max((n)/5,1));

    return 0;
}