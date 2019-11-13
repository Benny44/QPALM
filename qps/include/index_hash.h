#ifndef INDEX_HASH_H
#define INDEX_HASH_H

#include "global_opts.h"

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

struct index_table *create_index_table(c_int size);

c_int hashcode(struct index_table *t,char* key);

void insert(struct index_table *t, char* key, c_int index, char sign);

struct node* lookup(struct index_table *t,char* key);

void print_table(struct index_table *t, c_int size);

void free_index_table(struct index_table *t, c_int size);


#endif