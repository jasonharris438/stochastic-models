#include "stochastic_models/likelihood/general_linear_online.h"

#include <cmath>
#include <stdexcept>
GeneralLinearUpdater::GeneralLinearUpdater(const double mu_numerator,
                                           const double mu_denominator,
                                           const double sigma_state_value,
                                           const double last_observation,
                                           const double n_obs)
    : mu_numerator(mu_numerator),
      mu_denominator(mu_denominator),
      sigma_state_value(sigma_state_value),
      last_observation(last_observation),
      n_obs(n_obs) {}

const double GeneralLinearUpdater::initialSigmaState(const double& sigma,
                                                     const int& n_obs) {
    return std::pow(n_obs * sigma, 2);
}

void GeneralLinearUpdater::applyInternalMuUpdate(
    const double mu_numerator_new, const double mu_denominator_new,
    const double new_observation) {
    mu_numerator = mu_numerator_new;
    mu_denominator = mu_denominator_new;
    last_observation = new_observation;
}

void GeneralLinearUpdater::applyInternalSigmaUpdate(
    const double sigma_state_value_new, const int n_obs_new) {
    sigma_state_value = sigma_state_value_new;
    n_obs = n_obs_new;
}

const double GeneralLinearUpdater::updateMu(const double value) {
    const double mu_numerator_new = last_observation * value + mu_numerator;
    const double mu_denominator_new = std::pow(value, 2) + mu_denominator;

    double return_val{};
    if (mu_denominator_new == 0 || mu_numerator_new == 0) {
        return_val = 0;
    } else {
        return_val = std::log(mu_numerator_new / mu_denominator_new);
    }
    applyInternalMuUpdate(mu_numerator_new, mu_denominator_new, value);
    return return_val;
}

const double GeneralLinearUpdater::updateSigma(const double& value,
                                               const double& expected_x) {
    if (n_obs < 0) {
        throw std::invalid_argument("n_obs must be greater than 0.");
    }

    const double n_obs_new = n_obs + 1;
    const double n_obs_ratio = n_obs / n_obs_new;
    const double sigma_state_value_new =
        sigma_state_value + (n_obs_ratio * std::pow(value - expected_x, 2));
    const double return_val = std::sqrt(sigma_state_value_new / (n_obs_new));
    applyInternalSigmaUpdate(sigma_state_value_new, n_obs_new);
    return return_val;
}
