#ifndef SCALING_H
# define SCALING_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

// Functions to scale problem data
# include "types.h"
# include "lin_alg.h"
# include "constants.h"


/**
 * Scale problem matrices
 * @param  work Workspace
 * @return      exitflag
 */
void scale_data(QPALMWorkspace *work);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef SCALING_H