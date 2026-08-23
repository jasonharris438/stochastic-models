#ifndef STOCHASTIC_MODELS_NUMERIC_UTILS_INTEGRATION_H
#define STOCHASTIC_MODELS_NUMERIC_UTILS_INTEGRATION_H
#include "stochastic_models/numeric_utils/types.h"

#include <gsl/gsl_integration.h>

/**
 * @file
 * @brief RAII helpers and wrappers around GSL integration routines.
 */

/**
 * @brief RAII wrapper for a GSL integration workspace.
 *
 * This helper owns a gsl_integration_workspace pointer and ensures the
 * GSL workspace instance is correctly freed when the object is destroyed.
 */
class IntegrationState {
public:
  gsl_integration_workspace* workspace;
  /**
   * @brief Construct the RAII wrapper around an existing workspace reference.
   *
   * The constructor takes a reference to a GSL workspace that is managed by
   * this wrapper. The semantics mirror those in the implementation file.
   */
  IntegrationState(gsl_integration_workspace& w);
  /**
   * @brief Free the underlying GSL workspace if present.
   */
  ~IntegrationState();
};

/**
 * @brief Integrates fn over [lower, upper] with adaptive quadrature.
 *
 * The routine splits the interval where the error estimate is largest,
 * until the total estimate meets the tolerance.
 *
 * @param fn Function pointer conforming to ModelFunc that computes f(x).
 * @param model Opaque pointer passed to the function; used to carry model
 *              parameters or context.
 * @param lower Lower bound of integration (passed by reference to allow
 *              adaptive routines to modify it in some callers).
 * @param upper Upper bound of integration.
 * @return const double Value of the integral over [lower, upper].
 * @throws NonFiniteResultError If the integral value is NaN or Inf.
 * @throws IntegrationToleranceError If the error estimate is greater than
 * the accepted bound.
 */
const double
adaptiveIntegration(ModelFunc fn, void* model, double& lower, double& upper);
/**
 * @brief Integrates fn over [lower, +inf) with adaptive quadrature.
 *
 * The routine maps the semi-infinite interval onto a finite one, then
 * splits it where the error estimate is largest, until the total estimate
 * meets the tolerance.
 *
 * @param fn Function pointer conforming to ModelFunc that computes f(x).
 * @param model Opaque pointer passed to the function; used to carry model
 *              parameters or context.
 * @param lower Lower bound of the semi-infinite integral.
 * @return const double Value of the integral over [lower, +inf).
 * @throws NonFiniteResultError If the integral value is NaN or Inf.
 * @throws IntegrationToleranceError If the error estimate is greater than
 * the accepted bound.
 */
const double
semiInfiniteIntegrationUpper(ModelFunc fn, void* model, double& lower);
#endif // STOCHASTIC_MODELS_NUMERIC_UTILS_INTEGRATION_H
