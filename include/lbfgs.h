#ifndef LBFGS_H
# define LBFGS_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"

/**
 * Main lbfgs function
 * Sets work->d to the direction calculated by lbfgs
 * @param work Workspace
 */
void lbfgs_set_direction(QPALMWorkspace *work);

/**
 * Updates Sbuffer, Ybuffer, YSbuffer, currmem and curridx
 * @param work Workspace
 */
void lbfgs_update_buffers(QPALMWorkspace *work);

/**
 * Executes the two-loop lbfgs recursion. The result is stored in work->d
 * @param work Workspace
 */
void lbfgs_two_loop(QPALMWorkspace *work);


# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef SCALING_H