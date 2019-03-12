#ifndef VALIDATE_H
# define VALIDATE_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

# include "types.h"

/***********************************************************
* Validation of data and settings * *
***********************************************************/

/**
 * Validate problem data
 * @param  data QPALMData to be validated
 * @return      Exitflag to check
 */
c_int validate_data(const QPALMData *data);


/**
 * Validate problem settings
 * @param  settings QPALMSettings to be validated
 * @return      Exitflag to check
 */
c_int validate_settings(const QPALMSettings *settings);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef AUXIL_H