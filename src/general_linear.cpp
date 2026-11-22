#include "stochastic_models/sde/general_linear.h"

#include "stochastic_models/distributions/gaussian.h"

#include <cmath>
#include <limits>
#include <stdexcept>
GeneralLinearModel::GeneralLinearModel()
    : GeneralLinearModel::GeneralLinearModel(0.0, 1.0) {}
GeneralLinearModel::GeneralLinearModel(const double mu, const double sigma)
    : mu(mu), sigma(sigma) {
  dist = new GaussianDistribution();
}
GeneralLinearModel::GeneralLinearModel(
    const double mu, const double sigma, const std::uint64_t seed
)
    : mu(mu), sigma(sigma) {
  dist = new GaussianDistribution(0.0, 1.0, seed);
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
  if (mu >= 0) {
    throw std::domain_error(
        "The linear SDE has no stationary variance for mu >= 0."
    );
  }
  return std::pow(sigma, 2) / (-2 * mu);
}
const double GeneralLinearModel::getConditionalVariance(const double t) const {
  const double exponent = 2 * mu * t;
  // A subnormal exponent has lost precision, and the limit is exact there.
  if (mu == 0.0 || std::abs(exponent) < std::numeric_limits<double>::min()) {
    return std::pow(sigma, 2) * t;
  }
  // expm1 avoids catastrophic cancellation in exp(2 mu t) - 1 for small
  // mu * t.
  const double growth = std::expm1(exponent);
  // Dividing first keeps the product with sigma^2 from underflowing.
  const double quotient = growth / (2 * mu);
  if (std::isinf(quotient)) {
    // Multiplying first keeps the quotient from overflowing early.
    return std::pow(sigma, 2) * growth / (2 * mu);
  }
  return std::pow(sigma, 2) * quotient;
}
std::vector<double> GeneralLinearModel::Simulate(
    const double start, const unsigned int& size, const double t
) const {
  if (size < 1) {
    throw std::invalid_argument(
        "GeneralLinearModel::Simulate requires size >= 1."
    );
  }
  if (!std::isfinite(t) || t <= 0.0) {
    throw std::invalid_argument(
        "GeneralLinearModel::Simulate requires a finite t > 0."
    );
  }
  const std::vector<double> distribution_draws = (*dist).sample(size);
  std::vector<double> vec = {start};

  for (unsigned int n{}; n < size; n++) {
    const double sample = coreEquation(vec[n], distribution_draws[n], t);
    vec.push_back(sample);
  }

  return vec;
}
const double GeneralLinearModel::coreEquation(
    const double& x, const double& noise, const double t
) const {
  const double exp_mu_t = std::exp(mu * t);
  const double diffusion_sd = std::sqrt(getConditionalVariance(t));
  return x * exp_mu_t + diffusion_sd * noise;
}
