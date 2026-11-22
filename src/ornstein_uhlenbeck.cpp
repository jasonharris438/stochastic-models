#include "stochastic_models/sde/ornstein_uhlenbeck.h"

#include "stochastic_models/distributions/gaussian.h"

#include <cmath>
#include <limits>
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
 * @brief Seeded constructor for reproducible simulation.
 *
 */
OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel(
    const double mu,
    const double alpha,
    const double sigma,
    const std::uint64_t seed
)
    : mu(mu), alpha(alpha), sigma(sigma) {
  dist = new GaussianDistribution(0.0, 1.0, seed);
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
  const double exponent = 2 * alpha * t;
  // A subnormal exponent has lost precision, and the limit is exact there.
  if (alpha == 0.0 || std::abs(exponent) < std::numeric_limits<double>::min()) {
    return std::pow(sigma, 2) * t;
  }
  // expm1 avoids catastrophic cancellation in 1 - exp(-2 alpha t) for small
  // alpha * t.
  const double decay = -std::expm1(-exponent);
  // Dividing first keeps the product with sigma^2 from underflowing.
  const double quotient = decay / (2 * alpha);
  if (std::isinf(quotient)) {
    // Multiplying first keeps the quotient from overflowing early.
    return std::pow(sigma, 2) * decay / (2 * alpha);
  }
  return std::pow(sigma, 2) * quotient;
}
std::vector<double> OrnsteinUhlenbeckModel::Simulate(
    const double start, const unsigned int& size, const double t
) const {
  if (size < 1) {
    throw std::invalid_argument(
        "OrnsteinUhlenbeckModel::Simulate requires size >= 1."
    );
  }
  if (!std::isfinite(t) || t <= 0.0) {
    throw std::invalid_argument(
        "OrnsteinUhlenbeckModel::Simulate requires a finite t > 0."
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
    const double& x, const double& noise, const double t
) const {
  const double delta{std::exp(-alpha * t)};
  const double diffusion_sd = std::sqrt(getConditionalVariance(t));
  return (x * delta) + (mu * (1 - delta)) + (diffusion_sd * noise);
}
