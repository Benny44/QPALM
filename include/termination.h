#ifndef TERMINATION_H
# define TERMINATION_H

# ifdef __cplusplus
extern "C" {
# endif // ifdef __cplusplus

#include "types.h"

c_int check_termination(QPALMWorkspace *work);

void calculate_residuals_and_tolerances(QPALMWorkspace *work);

void calculate_primal_residual(QPALMWorkspace *work);

void calculate_dual_residual(QPALMWorkspace *work);

void calculate_primal_tolerance(QPALMWorkspace *work);

void calculate_dual_tolerances(QPALMWorkspace *work);

c_int is_solved(QPALMWorkspace *work);

c_int is_primal_infeasible(QPALMWorkspace *work);

c_int is_dual_infeasible(QPALMWorkspace *work);

c_int check_subproblem_termination(QPALMWorkspace *work);

void store_solution(QPALMWorkspace *work);



# ifdef __cplusplus
}
# endif // ifdef __cplusplus

#endif // ifndef TERMINATION_H