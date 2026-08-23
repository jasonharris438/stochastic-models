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
 * @brief Finds a root of fn inside [lower, upper] with Brent's method.
 *
 * Brent's method keeps an interval that brackets the root and shrinks it
 * each iteration. It stops when the interval is inside the tolerance.
 *
 * @param fn Function pointer to the scalar function whose root is sought.
 * @param model Opaque model/context pointer passed to fn. Contains model
 * instance that is being used.
 * @param lower Lower bound of the bracketing interval (may be updated).
 * @param upper Upper bound of the bracketing interval (may be updated).
 * @return const double The converged root value.
 * @throws std::invalid_argument If lower >= upper.
 * @throws NoMemoryError If solver allocation fails.
 * @throws RootNotBracketedError If fn does not change sign over the interval.
 * @throws NoSolutionError If fn cannot be evaluated at an iteration point.
 * @throws SolverConvergenceError If the iteration limit is reached before
 * the interval test passes.
 * @throws NonFiniteResultError If the converged value is NaN or Inf.
 */
const double
brentSolver(ModelFunc fn, void* model, double& lower, double& upper);
#endif // STOCHASTIC_MODELS_NUMERIC_UTILS_SOLVERS_H
