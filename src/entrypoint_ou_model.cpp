/**
 * @file entrypoint_ou_model.cpp
 * @brief Entry point module for Ornstein Uhlenbeck process calculation
 * functions.
 *
 */

#include "stochastic_models/entrypoints/ou_model.h"
#include "stochastic_models/hitting_times/hitting_time_density.h"
#include "stochastic_models/hitting_times/hitting_time_ornstein_uhlenbeck.h"
#include "stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h"
#include "stochastic_models/likelihood/ornstein_uhlenbeck_online.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"

const std::vector<double> simulateOrnsteinUhlenbeck(
    const double& mu,
    const double& alpha,
    const double& sigma,
    const double start,
    const unsigned int& size,
    const unsigned int& t
) {
  // Create the Ornstein-Uhlenbeck process model
  OrnsteinUhlenbeckModel model(mu, alpha, sigma);

  // Simulate the process
  return model.Simulate(start, size, t);
}

const double hittingTimeDensityOrnsteinUhlenbeck(
    double x,
    const double mu,
    const double alpha,
    const double sigma,
    double first,
    double second
) {
  // Create core model as pointer to void because that type is required by
  // GSL.
  void* model = new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);

  // The ModelFunc is provided to determine which function to provide to
  // GSL for integration.
  ModelFunc fn = &integrateHittingTimeDensity;

  // Calculate the hitting time density at point x given the arguments
  // first and second.
  const double value = hittingTimeDensity(x, fn, model, first, second);

  // Then need to cast model back to OrnsteinUhlenbeckModel pointer to
  // clean up.
  HittingTimeOrnsteinUhlenbeck* ptr =
      static_cast<HittingTimeOrnsteinUhlenbeck*>(model);
  delete ptr;
  ptr = nullptr;
  model = nullptr;
  return value;
}
const std::unordered_map<std::string, const double>
ornsteinUhlenbeckMaximumLikelihood(const std::vector<double> vec) {
  // Generate likelihood calculator and generate estimate of mu, alpha,
  // and sigma.
  OrnsteinUhlenbeckLikelihood likelihood = OrnsteinUhlenbeckLikelihood();
  const OrnsteinUhlenbeckLikelihoodComponents components =
      likelihood.calculateComponents(vec);
  const OrnsteinUhlenbeckParameters parameters =
      likelihood.calculateParameters(components);
  const std::unordered_map<std::string, const double> key_value_pairs{
      {"mu", parameters.mu},
      {"alpha", parameters.alpha},
      {"sigma", parameters.sigma}
  };
  return key_value_pairs;
}

const std::vector<double> updateOuModel(
    const double mu,
    const double alpha,
    const double sigma,
    const double lead_sum,
    const double lag_sum,
    const double lead_sum_squared,
    const double lag_sum_squared,
    const double lead_lag_sum_product,
    const uint32_t n_obs,
    const double new_observation,
    const double last_observation

) {
  OrnsteinUhlenbeckUpdater updater(
      OrnsteinUhlenbeckLikelihoodComponents{
          lead_sum, lag_sum, lead_sum_squared, lag_sum_squared,
          lead_lag_sum_product, n_obs
      },
      OrnsteinUhlenbeckParameters{mu, alpha, sigma}
  );

  const OrnsteinUhlenbeckParameters params =
      updater.updateState(new_observation, last_observation);

  return {params.mu, params.alpha, params.sigma};
}
