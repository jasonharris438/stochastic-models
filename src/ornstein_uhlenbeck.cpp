#include "stochastic_models/sde/ornstein_uhlenbeck.h"

#include "stochastic_models/distributions/gaussian.h"

#include <cmath>
#include <stdexcept>
/**
 * @brief No args constructor delegates to main constructor.
 *
 */
OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel()
    : OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel(0.0, 0.0, 1.0) {}
/**
 * @brief Main constructor.
 *
 */
OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel(
    const double mu, const double alpha, const double sigma
)
    : mu(mu), alpha(alpha), sigma(sigma) {
  dist = new GaussianDistribution();
}
/**
 * @brief Copy constructor.
 *
 * @param other The OrnsteinUhlenbeckModel to copy.
 */
OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel(
    const OrnsteinUhlenbeckModel& other
)
    : mu(other.mu), alpha(other.alpha), sigma(other.sigma) {
  dist = new GaussianDistribution();
}
/**
 * @brief Destructor.
 *
 */
OrnsteinUhlenbeckModel::~OrnsteinUhlenbeckModel() {
  delete dist;
  dist = nullptr;
}
const OrnsteinUhlenbeckModel* OrnsteinUhlenbeckModel::clone() const {
  return new OrnsteinUhlenbeckModel(*this);
}
const double OrnsteinUhlenbeckModel::getMean() const {
  return mu;
}
const double OrnsteinUhlenbeckModel::getUnconditionalVariance() const {
  return std::pow(sigma, 2) / (2 * alpha);
}
const double
OrnsteinUhlenbeckModel::getConditionalVariance(const double t) const {
  if (std::abs(alpha) < 1e-12) {
    return std::pow(sigma, 2) * t;
  }
  // expm1 avoids catastrophic cancellation in 1 - exp(-2 alpha t) for small
  // alpha * t; the guard above still handles the 0/0 case at alpha == 0.
  return std::pow(sigma, 2) * -std::expm1(-2 * alpha * t) / (2 * alpha);
}
std::vector<double> OrnsteinUhlenbeckModel::Simulate(
    const double start, const unsigned int& size, const unsigned int& t
) const {
  if (size < 1) {
    throw std::invalid_argument(
        "OrnsteinUhlenbeckModel::Simulate requires size >= 1."
    );
  }
  const std::vector<double> distribution_draws = (*dist).sample(size - 1);
  std::vector<double> vec = {start};

  for (double val : distribution_draws) {
    const double sample = coreEquation(vec.back(), val, t);
    vec.push_back(sample);
  }

  return vec;
}
const double OrnsteinUhlenbeckModel::coreEquation(
    const double& x, const double& noise, const unsigned int& t
) const {
  const double delta{std::exp(-alpha * t)};
  const double diffusion_sd = std::sqrt(getConditionalVariance(t));
  return (x * delta) + (mu * (1 - delta)) + (diffusion_sd * noise);
}
