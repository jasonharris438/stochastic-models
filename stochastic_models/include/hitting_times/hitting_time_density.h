#ifndef _HITTING_TIME_DENSITY_H_
#define _HITTING_TIME_DENSITY_H_
#include "numeric_utils/types.h"

/**
 * @brief Implements the function @f[ S(x) @f] that is used to produce the
 * hitting time density estimate.
 *
 * @param x The point at which to perform integration of the hitting time
 * density.
 * @param model A pointer to the model instance that is being used to
 * parameterize the integration.
 * @return double
 */
double integrate_hitting_time_density(double x, void* model);
/**
 * @brief Implements the function @f[ F(x) @f] that is used in the optimal mean
 * reversion trading problem.
 *
 * @param x The point at which to perform integration of @f[ F(x) @f].
 * @param model A pointer to the model instance that is being used to
 * parameterize the integration.
 * @return double
 */
double integrate_mean_reversion_f(double x, void* model);
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
const double hitting_time_density(double& x, model_func fn, void* model,
                                  double& first, double& second);
#endif  // _HITTING_TIME_DENSITY_H_
