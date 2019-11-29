#include "cholmod.h"
#include "qpalm.h"
#include "cholmod_interface.h"

cholmod_sparse *python_allocate_cholmod_sparse(size_t m, size_t n, size_t nzmax);
QPALMSettings *qpalm_malloc_settings(void);
QPALMData *qpalm_malloc_data(void);
void qpalm_free_settings(QPALMSettings *settings);
void qpalm_free_data(QPALMData *data);