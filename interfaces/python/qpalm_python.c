#include "qpalm.h"
#ifdef USE_LADEL
#include "ladel.h"
#elif defined USE_CHOLMOD
#include "cholmod.h"
#endif

solver_sparse *python_allocate_sparse(size_t m, size_t n, size_t nzmax) {
  solver_common common, *c;
  c = &common;
  solver_sparse *M;
  #ifdef USE_LADEL
  M = ladel_sparse_alloc(m, n, nzmax, UNSYMMETRIC, TRUE, FALSE);
  #elif defined USE_CHOLMOD
  CHOLMOD(start)(c);
  cholmod_set_settings(c);
  M = CHOLMOD(allocate_sparse)(m, n, nzmax, TRUE, TRUE, 0, CHOLMOD_REAL, c);
  CHOLMOD(finish)(c);
  #endif
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
    solver_common common, *c;
    c = &common;
    if (data) {
      #ifdef USE_LADEL
      data->Q = ladel_sparse_free(data->Q);
      data->A = ladel_sparse_free(data->A);
      #elif defined USE_CHOLMOD
      CHOLMOD(start)(c);
      if (data->Q) CHOLMOD(free_sparse)(&data->Q, c);
      if (data->A) CHOLMOD(free_sparse)(&data->A, c);
      CHOLMOD(finish)(c);
      #endif 
      
      if (data->q) c_free(data->q);

      if (data->bmin) c_free(data->bmin);

      if (data->bmax) c_free(data->bmax);
      c_free(data);
    }
}