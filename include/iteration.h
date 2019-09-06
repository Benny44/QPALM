/**
 * @file iteration.h
 * @author Ben Hermans
 * @brief QPALM main solver API.
 * @details This file contains the functions that make up the qpalm algorithm (the functions that qpalm_solve will use). 
 * These include the computation of the residuals at the start of the iteration, 
 * the update of the primal variables in an inner iteration, 
 * the update of the penalty factors and dual variables in an outer iteration, etc.
 */

#ifndef ITERATION_H
#define ITERATION_H

#include "types.h"
#include "global_opts.h"

void compute_residuals(QPALMWorkspace* work);

/**
 * Initialize penalty factors from initial x
 * 
 * The formula used here can be found in \cite birgin2014practical.  
 * @param work Workspace
 */
void initialize_sigma(QPALMWorkspace *work);

void update_sigma(QPALMWorkspace* work);

void update_gamma(QPALMWorkspace* work);

void update_primal_iterate(QPALMWorkspace *work);

c_float compute_objective(QPALMWorkspace *work);

#endif