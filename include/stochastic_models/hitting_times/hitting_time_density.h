#ifndef HITTING_TIME_DENSITY_H
#define HITTING_TIME_DENSITY_H
#include "stochastic_models/numeric_utils/types.h"

/**
 * @brief Implements the function @f[ S(x) @f] that is used to produce the
 * hitting time density estimate.
 *
 * @param x The point at which to perform integration of the hitting time
 * density.
 * @param model A pointer to the model instance that is being used to
 * parameterize the integration. This is statically cast to a
 * HittingTimeOrnsteinUhlenbeck instance.
 * @return double
 */
double integrateHittingTimeDensity(double x, void* model);
/**
 * @brief The hitting time density function @f[ \frac{S(x) - S(b)}{S(a) - S(b)}
 * @f].
 *
 * @param x Point at which to evaluate the hitting time density function.
 * @param fn A pointer to the function to integrate.
 * @param model A pointer to the model instance that is being used
 * to perform integration.
 * @param first The first hitting time.
 * @param second The second hitting time.
 * @return const double. The value of the hitting time density function at @f[
 * x@f].
 */
const double hittingTimeDensity(double& x, ModelFunc fn, void* model,
                                double& first, double& second);
#endif  // HITTING_TIME_DENSITY_H
