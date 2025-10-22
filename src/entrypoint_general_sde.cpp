/**
 * @file entrypoint_general_sde.cpp
 * @brief Entry point module for General Linear SDE model calculation
 * functions.
 *
 */

#include "stochastic_models/entrypoints/general_sde.h"
#include "stochastic_models/likelihood/general_linear_online.h"

#include <cmath>

const std::vector<double> updateGeneralLinearSDEModel(
    const double mu,
    const double sigma,
    const double new_observation,
    const double last_observation,
    const double lag_squared,
    const double lead_lag_inner_product,
    const double squared_error,
    const uint32_t n_obs
) {
  GeneralLinearUpdater updater(
      GeneralLinearLikelihoodComponents{
          lag_squared, lead_lag_inner_product, squared_error, n_obs
      },
      GeneralLinearParameters{mu, sigma}
  );

  const GeneralLinearParameters params =
      updater.updateState(new_observation, last_observation);
  return {params.mu, params.sigma};
}

