#ifndef QPALM_H
#define QPALM_H

#include "types.h"

/********************
* Main Solver API  *
********************/

/**
 * @name Main solver API
 * @{
 */


/**
 * Set default settings from constants.h file
 * assumes settings already allocated in memory
 * @param settings settings structure
 */
void qpalm_set_default_settings(QPALMSettings *settings);


/**
 * Initialize QPALM solver allocating memory.
 *
 * All the inputs must be already allocated in memory before calling.
 *
 * It performs:
 * - data and settings validation
 * - problem data scaling
 *
 * NB: This is the only function that allocates dynamic memory and is not used
 *during code generation
 *
 * @param  data         Problem data
 * @param  settings     Solver settings
 * @return              Solver environment
 */
QPALMWorkspace* qpalm_setup(const QPALMData *data,
                          QPALMSettings   *settings);



/**
 * Solve quadratic program
 *
 * The final solver information is stored in the \a work->info  structure
 *
 * The solution is stored in the  \a work->solution  structure
 *
 * If the problem is primal infeasible, the certificate is stored
 * in \a work->delta_y
 *
 * If the problem is dual infeasible, the certificate is stored in \a
 * work->delta_x
 *
 * @param  work Workspace allocated
 */
void qpalm_solve(QPALMWorkspace *work);



/**
 * Cleanup workspace by deallocating memory
 *
 * This function is not used in code generation
 * @param  work Workspace
 */
void qpalm_cleanup(QPALMWorkspace *work);

#endif