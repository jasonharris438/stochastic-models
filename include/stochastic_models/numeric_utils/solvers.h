#ifndef SOLVER_H
#define SOLVER_H
#include <gsl/gsl_roots.h>

#include "stochastic_models/numeric_utils/types.h"

/**
 * @brief RAII wrapper for GSL root solver.
 *
 * This class manages the lifecycle of a GSL root solver
 * by allocating and deallocating the solver as needed. It ensures
 * that memory is properly freed when the object goes out of scope,
 */
class BrentSolverState {
   public:
    gsl_root_fsolver* fsolver;
    const gsl_root_fsolver_type* fsolver_type;
    /**
     * @brief Constructor that allocates memory for the Brent root solver.
     * It initializes the solver type and allocates the solver workspace.
     */
    BrentSolverState();
    /**
     * @brief Destructor that frees the memory allocated for the Brent
     * root solver. It checks if the solver pointer is not null before
     * freeing it using native GSL functions.
     */
    ~BrentSolverState();
};

/**
 * @brief Uses the Brent optimization method to solve for the root of a
 * function.
 *
 * @param fn A pointer to the function to solve.
 * @param model A pointer to the model instance that is being used
 * to perform root finding.
 * @param lower A reference to the lower bound of the root finding interval.
 * @param upper A reference to the upper bound of the root finding interval.
 * @return const double. The value of the root of the function.
 */
const double brentSolver(ModelFunc fn, void* model, double& lower,
                         double& upper);
#endif  // SOLVER_H
