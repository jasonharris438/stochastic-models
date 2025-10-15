#include "stochastic_models/likelihood/general_linear_online.h"

#include "stochastic_models/exceptions/errors.h"

#include <cmath>
GeneralLinearUpdater::GeneralLinearUpdater(
    const double mu_numerator,
    const double mu_denominator,
    const double sigma_state_value,
    const double last_observation,
    const int n_obs
)
    : mu_numerator(mu_numerator), mu_denominator(mu_denominator),
      sigma_state_value(sigma_state_value), last_observation(last_observation),
      n_obs(n_obs) {}

const double
GeneralLinearUpdater::getExpectedValue(const double& observation) const {
  return std::exp(getMu()) * observation;
}

const double GeneralLinearUpdater::getMu() const {
  return std::log(mu_numerator / mu_denominator);
}
const double GeneralLinearUpdater::getSigma() const {
  if (n_obs <= 1) {
    throw InvalidNumberObservationsError("n_obs must be greater than 1.");
  }
  return std::sqrt(sigma_state_value / n_obs);
}

void GeneralLinearUpdater::applyInternalMuUpdate(
    const double mu_numerator_new, const double mu_denominator_new
) {
  mu_numerator = mu_numerator_new;
  mu_denominator = mu_denominator_new;
}

void GeneralLinearUpdater::applyInternalSigmaUpdate(
    const double sigma_state_value_new
) {
  sigma_state_value = sigma_state_value_new;
}
void GeneralLinearUpdater::internalStateUpdate(
    const int n_obs_new, const double new_observation
) {
  n_obs = n_obs_new;
  last_observation = new_observation;
}

const double GeneralLinearUpdater::updateMu(const double& value) {
  const double mu_numerator_new = last_observation * value + mu_numerator;
  const double mu_denominator_new =
      std::pow(last_observation, 2) + mu_denominator;

  double return_val{};
  if (mu_denominator_new == 0 || mu_numerator_new == 0) {
    return_val = 0;
  } else {
    return_val = std::log(mu_numerator_new / mu_denominator_new);
  }
  applyInternalMuUpdate(mu_numerator_new, mu_denominator_new);
  return return_val;
}

const double GeneralLinearUpdater::updateSigma(const double& value) {
  if (n_obs < 0) {
    throw InvalidNumberObservationsError("n_obs must be greater than 0.");
  }

  const double n_obs_ratio = n_obs / (n_obs + 1);
  const double error = value - getExpectedValue(last_observation);
  const double sigma_state_value_new =
      sigma_state_value + (n_obs_ratio * std::pow(error, 2));
  const double return_val = std::sqrt(sigma_state_value_new / (n_obs + 1));
  applyInternalSigmaUpdate(sigma_state_value_new);
  return return_val;
}
