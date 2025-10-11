#include "stochastic_models/sde/general_linear.h"

#include <cmath>
GeneralLinearModel::GeneralLinearModel()
    : GeneralLinearModel::GeneralLinearModel(0.0, 1.0) {}
GeneralLinearModel::GeneralLinearModel(const double mu, const double sigma)
    : mu(mu), sigma(sigma) {
  dist = new GaussianDistribution();
}
GeneralLinearModel::GeneralLinearModel(const GeneralLinearModel& other)
    : mu(other.mu), sigma(other.sigma) {
  dist = new GaussianDistribution();
}
const GeneralLinearModel* GeneralLinearModel::clone() const {
  return new GeneralLinearModel(*this);
}
GeneralLinearModel::~GeneralLinearModel() {
  delete dist;
  dist = nullptr;
}
const double GeneralLinearModel::getMean() const {
  return 0.0;
}
const double GeneralLinearModel::getUnconditionalVariance() const {
  if (mu == 0) {
    return 0;
  }
  return (std::pow(sigma, 2) / (2 * mu)) * (std::exp(2 * mu) - 1);
}
const double GeneralLinearModel::getConditionalVariance() const {
  if (mu == 0) {
    return 0;
  }
  return ((2 * sigma * mu) / (std::exp(2 * mu) - std::exp(mu)));
}
std::vector<double> GeneralLinearModel::Simulate(
    const double start, const unsigned int& size, const unsigned int& t
) const {
  const std::vector<double> distribution_draws = (*dist).sample(size);
  std::vector<double> vec = {start};

  for (unsigned int n{}; n < size; n++) {
    const double sample = coreEquation(vec[n], distribution_draws[n], t);
    vec.push_back(sample);
  }

  return vec;
}
const double GeneralLinearModel::coreEquation(
    const double& x, const double& noise, const unsigned int& t
) const {
  const double exp_mu_t = std::exp(mu * t);
  return x * exp_mu_t + exp_mu_t * std::exp(-mu * t) * sigma * noise;
}
