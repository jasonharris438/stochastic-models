#ifndef STOCHASTIC_MODELS_HITTING_TIMES_HITTING_TIME_DENSITY_H
#define STOCHASTIC_MODELS_HITTING_TIMES_HITTING_TIME_DENSITY_H
#include "stochastic_models/numeric_utils/types.h"

/**
 * @file
 * @brief Helpers for computing hitting-time densities for OU models used by
 * the trading-level optimizers.
 */

/**
 * @brief Integrand function S(x) used by the hitting time density integrators.
 *
 * @param x Integration point.
 * @param model Opaque model pointer (expected to be
 * HittingTimeOrnsteinUhlenbeck*).
 * @return double Integrand value.
 */
double integrateHittingTimeDensity(double x, void* model);

/**
 * @brief Computes the normalized hitting-time density.
 *
 * The density is the integral of the core over [second, x] divided by the
 * integral over [second, first]. Exceptions from the integration wrappers
 * propagate to the caller.
 *
 * @param x Point at which to evaluate (may be modified by integrator routines).
 * @param fn Integrand function pointer.
 * @param model Opaque model pointer.
 * @param first Left boundary value.
 * @param second Right boundary value.
 * @return const double Evaluated hitting-time density.
 * @throws ZeroDivError If the denominator integral is zero.
 * @throws NonFiniteResultError If the density value is NaN or Inf.
 */
const double hittingTimeDensity(
    double& x, ModelFunc fn, void* model, double& first, double& second
);
#endif // STOCHASTIC_MODELS_HITTING_TIMES_HITTING_TIME_DENSITY_H
