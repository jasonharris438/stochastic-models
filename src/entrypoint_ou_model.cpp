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
  OrnsteinUhlenbeckLikelihood* likelihood = new OrnsteinUhlenbeckLikelihood();
  const std::unordered_map<std::string, const double> params{
      likelihood->calculate(vec)
  };
  delete likelihood;
  return params;
}
