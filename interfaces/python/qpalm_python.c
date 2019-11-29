#include "qpalm.h"
#include "global_opts.h"
#include "cholmod.h"
#include "constants.h"
#include "cholmod_interface.h"

cholmod_sparse *python_allocate_cholmod_sparse(size_t m, size_t n, size_t nzmax) {
  cholmod_common common, *c;
  c = &common;
  CHOLMOD(start)(c);
  cholmod_set_settings(c);
  cholmod_sparse *M = CHOLMOD(allocate_sparse)(m, n, nzmax, TRUE, TRUE, 0, CHOLMOD_REAL, c);
  CHOLMOD(finish)(c);
  return M;
}

QPALMSettings *python_allocate_settings(void) {
  return (QPALMSettings *) c_malloc(sizeof(QPALMSettings));
}

QPALMData *python_allocate_data(void) {
  return (QPALMData *) c_malloc(sizeof(QPALMData));
}

void python_free_settings(QPALMSettings *settings) {
  if (settings) c_free(settings);
}

void python_free_data(QPALMData *data) {
    cholmod_common common, *c;
    c = &common;
    if (data) {
      CHOLMOD(start)(c);

      if (data->Q) CHOLMOD(free_sparse)(&data->Q, c);

      if (data->A) CHOLMOD(free_sparse)(&data->A, c);

      CHOLMOD(finish)(c);

      if (data->q) c_free(data->q);

      if (data->bmin) c_free(data->bmin);

      if (data->bmax) c_free(data->bmax);
      c_free(data);
    }
}