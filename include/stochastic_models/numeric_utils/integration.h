#ifndef INTEGRATION_H
#define INTEGRATION_H
#include <gsl/gsl_integration.h>

#include "stochastic_models/numeric_utils/types.h"

/**
 * @brief RAII wrapper for GSL integration workspace.
 *
 * This class manages the lifecycle of a GSL integration workspace
 * by allocating and deallocating the workspace as needed. It ensures
 * that memory is properly freed when the object goes out of scope,
 */
class IntegrationState {
   public:
    gsl_integration_workspace* workspace;
    IntegrationState(gsl_integration_workspace& w);
    ~IntegrationState();
};

/**
 * @brief Integrates the function f over a given interval.
 *
 * @param fn A pointer to the function to integrate.
 * @param model A pointer to the model instance that is being used
 * to perform integration.
 * @param lower A reference to the lower bound of the integration interval.
 * @param upper A reference to the upper bound of the integration interval.
 * @return const double. The value of the integral of F.
 */
const double adaptiveIntegration(ModelFunc fn, void* model, double& lower,
                                 double& upper);
/**
 * @brief Integrates the function f over a semi-infinite interval [lower, +inf).
 *
 * @param fn A pointer to the function to integrate.
 * @param model A pointer to the model instance that is being used
 * to perform integration.
 * @param lower A reference to the lower bound of the integration interval.
 * @return const double. The value of the integral of F.
 */
const double semiInfiniteIntegrationUpper(ModelFunc fn, void* model,
                                          double& lower);
#endif  // INTEGRATION_H
