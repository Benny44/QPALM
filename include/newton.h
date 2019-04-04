#ifndef NEWTON_H
#define NEWTON_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"


/**
 * Main newton function
 * Sets work->d to the direction calculated by the semismooth Newton method
 * @param work Workspace
 */
void newton_set_direction(QPALMWorkspace *work);


/**
 * Computes and stores the set of active constraints
 * @param work Workspace
 */
void set_active_constraints(QPALMWorkspace *work);

/**
 * Determines the entering and leaving constraints
 * @param work Workspace
 */
void set_entering_leaving_constraints(QPALMWorkspace *work);

# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif //ifndef NEWTON_H