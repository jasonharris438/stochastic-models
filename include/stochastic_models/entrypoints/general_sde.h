#ifndef STOCHASTIC_MODELS_ENTRYPOINTS_GENERAL_SDE_H
#define STOCHASTIC_MODELS_ENTRYPOINTS_GENERAL_SDE_H

#include <cstdint>
#include <vector>

/**
 * @brief Handles updating the parameters of a General Linear SDE model given a
 * new observation and the current state of the model.
 *
 * The function arguments are granular to ensure efficiency when updating the
 * model state. These values should be managed by the calling code to maintain
 * consistency.
 *
 * @param observation The new observation to update the model with.
 * @param mu_numerator The current numerator value for mu calculation.
 * @param mu_denominator The current denominator value for mu calculation.
 * @param sigma_state_value The current state value for sigma calculation.
 * @param last_observation The last observation used in the model.
 * @param n_obs The current number of observations used in the model.
 * @return const std::vector<double> A vector containing the updated mu and
 * sigma values.
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
