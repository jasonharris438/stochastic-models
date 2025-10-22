#ifndef STOCHASTIC_MODELS_ENTRYPOINTS_GENERAL_SDE_H
#define STOCHASTIC_MODELS_ENTRYPOINTS_GENERAL_SDE_H

#include <cstdint>
#include <vector>

/**
 * @file
 * @brief Entrypoint used by the command-line example to perform a
 * single online update of the General Linear SDE model estimates.
 */

/**
 * @brief Update the model parameters and tracking statistics when a single new
 * observation is received.
 *
 * The function takes a minimal set of numeric accumulators in addition to the
 * current parameter estimates and returns a vector with updated mu and sigma
 * values (in that order).
 *
 * @param mu Current mu parameter estimate.
 * @param sigma Current sigma estimate.
 * @param new_observation New observed value X_t.
 * @param last_observation Previous observed value X_{t-1}.
 * @param lag_squared Current accumulated sum of lag^2.
 * @param lead_lag_inner_product Current accumulated lead-lag inner product.
 * @param squared_error Current accumulated squared error.
 * @param n_obs Current number of observations (before this update).
 * @return const std::vector<double> Vector of updated parameters: {mu, sigma}.
 */
const std::vector<double> updateGeneralLinearSDEModel(
    const double mu,
    const double sigma,
    const double new_observation,
    const double last_observation,
    const double lag_squared,
    const double lead_lag_inner_product,
    const double squared_error,
    const uint32_t n_obs
);

#endif // STOCHASTIC_MODELS_ENTRYPOINTS_GENERAL_SDE_H
