/**
 * @file entrypoint_general_sde.cpp
 * @brief Entry point module for General Linear SDE model calculation
 * functions.
 *
 */

#include "stochastic_models/entrypoints/general_sde.h"
#include "stochastic_models/likelihood/general_linear_online.h"

#include <cmath>

const std::vector<double> updateGeneralLinearSDE(
    const double observation,
    const double mu_numerator,
    const double mu_denominator,
    const double sigma_state_value,
    const double last_observation,
    const int n_obs
) {
  GeneralLinearUpdater updater(
      mu_numerator, mu_denominator, sigma_state_value, last_observation, n_obs
  );

  const double mu = updater.updateMu(observation);

  return {mu, updater.updateSigma(observation)};
}

