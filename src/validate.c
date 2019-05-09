/**
 * @file validate.c
 * @author Ben Hermans
 * @brief Validation of the user provided settings and data.
 * @details The assumptions on the settings can be found in the 
 * details of the QPALMSettings tab in data structures.
 */

#include "validate.h"
#include "lin_alg.h"
#include "constants.h"


/***********************************************************
* Validation of data and settings * *
***********************************************************/

c_int validate_data(const QPALMData *data) {
  

  if (!data) {
# ifdef PRINTING
    c_eprint("Missing data");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  // Lower and upper bounds
  size_t j;
  for (j = 0; j < data->m; j++) {
    if (data->bmin[j] > data->bmax[j]) {
# ifdef PRINTING
      c_eprint("Lower bound at index %d is greater than upper bound: %.4e > %.4e",
               (int)j, data->bmin[j], data->bmax[j]);
# endif /* ifdef PRINTING */
      return FALSE;
    }
  }
  return TRUE;
}


c_int validate_settings(const QPALMSettings *settings) {
  
  if (!settings) {
# ifdef PRINTING
    c_eprint("Missing settings!");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->max_iter <= 0) {
# ifdef PRINTING
    c_eprint("max_iter must be positive");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->eps_abs < 0) {
# ifdef PRINTING
    c_eprint("eps_abs must be nonnegative");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->eps_rel < 0) {
# ifdef PRINTING
    c_eprint("eps_rel must be nonnegative");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if ((settings->eps_rel == 0) && (settings->eps_abs == 0)) {
# ifdef PRINTING
    c_eprint("at least one of eps_abs and eps_rel must be positive");
# endif /* ifdef PRINTING */
    return FALSE;
  }

    if (settings->eps_abs_in < 0) {
# ifdef PRINTING
    c_eprint("eps_abs_in must be nonnegative");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->eps_rel_in < 0) {
# ifdef PRINTING
    c_eprint("eps_rel_in must be nonnegative");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if ((settings->eps_rel_in == 0) && (settings->eps_abs_in == 0)) {
# ifdef PRINTING
    c_eprint("at least one of eps_abs_in and eps_rel_in must be positive");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->rho <= 0 || settings->rho >= 1) {
# ifdef PRINTING
    c_eprint("rho must be positive and smaller than 1");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->eps_prim_inf < 0) {
# ifdef PRINTING
    c_eprint("eps_prim_inf must be nonnegative");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->eps_dual_inf < 0) {
# ifdef PRINTING
    c_eprint("eps_dual_inf must be nonnegative");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->theta > 1) {
# ifdef PRINTING
    c_eprint("theta must be smaller than ot equal 1");
# endif /* ifdef PRINTING */
    return FALSE;
  }

    if (settings->delta <= 1) {
# ifdef PRINTING
    c_eprint("delta must be greater than 1");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if (settings->tau_init <= 0) {
# ifdef PRINTING
    c_eprint("tau_init must be greater than 0");
# endif /* ifdef PRINTING */
    return FALSE;
  }

   if ((settings->proximal != 0) && (settings->proximal != 1)) {
# ifdef PRINTING
    c_eprint("proximal must be either 0 or 1");
# endif /* ifdef PRINTING */
    return FALSE;
  }

   if (settings->gamma_init <= 0) {
# ifdef PRINTING
    c_eprint("gamma_init must be positive");
# endif /* ifdef PRINTING */
    return FALSE;
  }

   if (settings->gamma_upd < 1) {
# ifdef PRINTING
    c_eprint("gamma update factor must be greater than or equal to 1");
# endif /* ifdef PRINTING */
    return FALSE;
  }

   if (settings->gamma_max < settings->gamma_init) {
# ifdef PRINTING
    c_eprint("gamma max must be greater than or equal to gamma");
# endif /* ifdef PRINTING */
    return FALSE;
  }

   if (settings->scaling < 0) {
# ifdef PRINTING
    c_eprint("scaling must be greater than or equal to zero");
# endif /* ifdef PRINTING */
    return FALSE;
  }

  if ((settings->warm_start != 0) && (settings->warm_start != 1)) {
# ifdef PRINTING
    c_eprint("warm_start must be either 0 or 1");
# endif /* ifdef PRINTING */
    return FALSE;
  }


  if ((settings->verbose != 0) && (settings->verbose != 1)) {
# ifdef PRINTING
    c_eprint("verbose must be either 0 or 1");
# endif /* ifdef PRINTING */
    return FALSE;
  }


  return TRUE;
}