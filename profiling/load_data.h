#include "types.h"
#include "cholmod.h"

void load_data(const char* name, QPALMData* data);
c_float *load_dense(FILE *fp, size_t n);
