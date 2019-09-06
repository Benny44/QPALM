/**
 * @file util.c
 * @author Ben Hermans
 * @brief Utility functions.
 * @details This file contains some utility functions, to copy the settings, initialize the penalty parameters,
 * initialize the iterates, update the solver status and time the algorithm.
 */

#include "util.h"
#include "lin_alg.h"
#include "global_opts.h"
/**********************
* Utility Functions  *
**********************/

void c_strcpy(char dest[], const char source[]) {
    size_t i;
    for(i = 0; (dest[i] = source[i]) != '\0'; i++);
}


QPALMSettings* copy_settings(const QPALMSettings *settings) {
    QPALMSettings *new = c_malloc(sizeof(QPALMSettings));

    // Copy settings
    new->max_iter     = settings->max_iter;     
    new->eps_abs      = settings->eps_abs;       
    new->eps_rel      = settings->eps_rel;       
    new->eps_abs_in   = settings->eps_abs_in;    
    new->eps_rel_in   = settings->eps_rel_in;    
    new->rho          = settings->rho;           
    new->eps_prim_inf = settings->eps_prim_inf;  
    new->eps_dual_inf = settings->eps_dual_inf; 
    new->theta        = settings->theta;         
    new->delta        = settings->delta;
    new->tau_init     = settings->tau_init;         
    new->proximal     = settings->proximal;       
    new->gamma_init   = settings->gamma_init;         
    new->gamma_upd    = settings->gamma_upd;     
    new->gamma_max    = settings->gamma_max;     
    new->scaling      = settings->scaling;    
    new->nonconvex    = settings->nonconvex;  
    new->verbose      = settings->verbose; 
    new->warm_start   = settings->warm_start;      

    return new;
}

void update_status(QPALMInfo *info, c_int status_val) {
    // Update status value
    info->status_val = status_val;

    // Update status string depending on status val
    switch (status_val)
    {
    case QPALM_SOLVED:
      c_strcpy(info->status, "solved");
      break;
    case QPALM_PRIMAL_INFEASIBLE:
      c_strcpy(info->status, "primal infeasible");
      break;
    case QPALM_DUAL_INFEASIBLE:
      c_strcpy(info->status, "dual infeasible");
      break;
    case QPALM_MAX_ITER_REACHED:
      c_strcpy(info->status, "maximum iterations reached");
      break;
    case QPALM_UNSOLVED:
      c_strcpy(info->status, "unsolved");
      break;
    case QPALM_ERROR:
      c_strcpy(info->status, "error");
      break;
    default:
      #ifdef PRINTING
        c_eprint("Unrecognised status value %d", status_val);
      #endif
      break;
    }
}

/*******************
* Timer Functions *
*******************/

#ifdef PROFILING

// Windows
# ifdef _WIN32

void qpalm_tic(QPALMTimer *t)
{
  QueryPerformanceFrequency(&t->freq);
  QueryPerformanceCounter(&t->tic);
}

c_float qpalm_toc(QPALMTimer *t)
{
  QueryPerformanceCounter(&t->toc);
  return (t->toc.QuadPart - t->tic.QuadPart) / (c_float)t->freq.QuadPart;
}

// Mac
# elif defined __APPLE__

void qpalm_tic(QPALMTimer *t)
{
  /* read current clock cycles */
  t->tic = mach_absolute_time();
}

c_float qpalm_toc(QPALMTimer *t)
{
  uint64_t duration; /* elapsed time in clock cycles*/

  t->toc   = mach_absolute_time();
  duration = t->toc - t->tic;

  /*conversion from clock cycles to nanoseconds*/
  mach_timebase_info(&(t->tinfo));
  duration *= t->tinfo.numer;
  duration /= t->tinfo.denom;

  return (c_float)duration / 1e9;
}

// Mac
# elif defined __MACH__

void qpalm_tic(QPALMTimer *t)
{
  /* read current clock cycles */
  t->tic = mach_absolute_time();
}

c_float qpalm_toc(QPALMTimer *t)
{
  uint64_t duration; /* elapsed time in clock cycles*/

  t->toc   = mach_absolute_time();
  duration = t->toc - t->tic;

  /*conversion from clock cycles to nanoseconds*/
  mach_timebase_info(&(t->tinfo));
  duration *= t->tinfo.numer;
  duration /= t->tinfo.denom;

  return (c_float)duration / 1e9;
}

// Linux
# elif defined __linux__  /* ifdef _WIN32 */
/* read current time */

void qpalm_tic(QPALMTimer *t)
{
  clock_gettime(CLOCK_MONOTONIC, &t->tic);
}

/* return time passed since last call to tic on this timer */
c_float qpalm_toc(QPALMTimer *t)
{
  struct timespec temp;

  clock_gettime(CLOCK_MONOTONIC, &t->toc);

  if ((t->toc.tv_nsec - t->tic.tv_nsec) < 0) {
    temp.tv_sec  = t->toc.tv_sec - t->tic.tv_sec - 1;
    temp.tv_nsec = 1e9 + t->toc.tv_nsec - t->tic.tv_nsec;
  } else {
    temp.tv_sec  = t->toc.tv_sec - t->tic.tv_sec;
    temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
  }
  return (c_float)temp.tv_sec + (c_float)temp.tv_nsec / 1e9;
}

# endif /* ifdef IS_WINDOWS */

#endif // If Profiling end
