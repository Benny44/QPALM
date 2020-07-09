#include "qpalm.h"

solver_sparse *python_allocate_sparse(size_t m, size_t n, size_t nzmax);

QPALMSettings *python_allocate_settings(void);

QPALMData *python_allocate_data(void);

void python_free_settings(QPALMSettings *settings);

void python_free_data(QPALMData *data);