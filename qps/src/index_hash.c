#include "index_hash.h"

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
        }
    }
    c_free(t->list);
    c_free(t);
}