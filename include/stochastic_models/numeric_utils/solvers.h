#ifndef STOCHASTIC_MODELS_NUMERIC_UTILS_SOLVERS_H
#define STOCHASTIC_MODELS_NUMERIC_UTILS_SOLVERS_H
#include "stochastic_models/numeric_utils/types.h"

#include <gsl/gsl_roots.h>

/**
 * @file
 * @brief Wrappers around GSL root-finding helpers.
 */

/**
 * @brief RAII wrapper for a GSL root solver (Brent method).
 *
 * Owns solver memory and ensures proper cleanup on destruction.
 */
class BrentSolverState {
public:
  gsl_root_fsolver* fsolver;
  const gsl_root_fsolver_type* fsolver_type;
  BrentSolverState();
  ~BrentSolverState();
};

/**
 * @brief Uses the Brent method to find a root of fn in [lower, upper].
 *
 * @param fn Function pointer to the scalar function whose root is sought.
 * @param model Opaque model/context pointer passed to fn. Contains model
 * instance that is being used.
 * @param lower Lower bound of the bracketing interval (may be updated).
 * @param upper Upper bound of the bracketing interval (may be updated).
 * @return const double Approximated root value.
 */
const double
brentSolver(ModelFunc fn, void* model, double& lower, double& upper);
#endif // STOCHASTIC_MODELS_NUMERIC_UTILS_SOLVERS_H
