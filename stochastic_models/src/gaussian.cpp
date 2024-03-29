#include "distributions/gaussian.h"

#include <time.h>

#include <cmath>
#include <random>
GaussianDistribution::GaussianDistribution(const double mu, const double sigma)
    : mu(mu), sigma(sigma) {}
GaussianDistribution::GaussianDistribution() : GaussianDistribution(0, 1.0) {}
const double GaussianDistribution::getMean() const {  // Returns class mean.
    return mu;
}
const double GaussianDistribution::getStd() const {  // Returns class std.
    return sigma;
}
double GaussianDistribution::erfGaussianCdf(
    const float x) {  // Core calculation for cdf evaluated at x.
    return (1.0 / 2.0) * (1.0 + std::erf(x));
}
const double GaussianDistribution::Cdf(
    const double& x) {  // Produces cdf for a given x.
    const double val = (x - mu) / (sigma * std::sqrt(2));
    return GaussianDistribution::erfGaussianCdf(val);
}
std::vector<double> GaussianDistribution::Sample(
    std::size_t size) const {  // Draws random samples from distribution.
    // Create uniform distribution instance to generate seed for normal
    // distribution
    std::uniform_int_distribution<> distrib(1, 99999999);

    // draw seed from uniform distribution
    const long int time_seed{static_cast<long int>(time(NULL))};
    std::mt19937 unif_gen(time_seed);
    int seed{distrib(unif_gen)};

    // instantiate normal distribution generator
    std::normal_distribution<> norm(mu, sigma);
    std::mt19937 norm_gen(seed);

    // generate size many samples
    std::vector<double> sample{};
    for (unsigned int i{0}; i < size; i++) {
        sample.push_back(norm(norm_gen));
    }

    return sample;
}
