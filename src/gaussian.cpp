#include "stochastic_models/distributions/gaussian.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace {
  // Mixes two random_device words into a 64-bit nondeterministic seed.
  std::uint64_t nondeterministic_seed() {
    std::random_device rd;
    return (static_cast<std::uint64_t>(rd()) << 32) ^
           static_cast<std::uint64_t>(rd());
  }
} // namespace
GaussianDistribution::~GaussianDistribution() {}
GaussianDistribution::GaussianDistribution(
    const double mu, const double sigma, const std::uint64_t seed
)
    : mu(mu), sigma(sigma), gen_(seed) {}
GaussianDistribution::GaussianDistribution(const double mu, const double sigma)
    : GaussianDistribution(mu, sigma, nondeterministic_seed()) {}
GaussianDistribution::GaussianDistribution() : GaussianDistribution(0.0, 1.0) {}
const double GaussianDistribution::getMean() const { // Returns class mean.
  return mu;
}
const double GaussianDistribution::getStd() const { // Returns class std.
  return sigma;
}
double GaussianDistribution::erfGaussianCdf(
    const float x
) { // Core calculation for cdf evaluated at x.
  return (1.0 / 2.0) * (1.0 + std::erf(x));
}
const double
GaussianDistribution::cdf(const double& x) { // Produces cdf for a given x.
  const double numerator = (x - mu);
  const double denominator = (sigma * std::sqrt(2));

  double val{};
  if (denominator == 0) {
    val = 0;
  } else {
    val = numerator / denominator;
  }
  return GaussianDistribution::erfGaussianCdf(val);
}
std::vector<double> GaussianDistribution::sample(
    const std::size_t& size
) const { // Draws random samples from the distribution.
  std::normal_distribution<double> norm(mu, sigma);
  std::vector<double> out(size);
  std::generate(out.begin(), out.end(), [&] { return norm(gen_); });
  return out;
}
