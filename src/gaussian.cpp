#include "gaussian.h"

#include <cmath>
#include <random>

GaussianDistribution::GaussianDistribution(const double mu, const double sigma) :
    mu(mu), sigma(sigma) { }

GaussianDistribution::GaussianDistribution() :
    GaussianDistribution(0, 1.0) { }

const double &GaussianDistribution::getMean() const {    // Returns const reference to class mean.
    const double& mu_ref = mu;
    return mu_ref;
}

const double &GaussianDistribution::getStd() const {    // Returns const reference to class std.
    const double& sigma_ref = sigma;
    return sigma_ref;
}

double GaussianDistribution::erfGaussianCdf(const float x){    // Core calculation for cdf evaluated at x.
    return (1.0/2.0) * (1.0 + std::erf(x));
}

double GaussianDistribution::Cdf(const double x){    // Produces cdf for a given x.
    const double val = (x - mu) / (sigma * std::sqrt(2));
    return GaussianDistribution::erfGaussianCdf(val);
}

double GaussianDistribution::Sample() const {   // Draws random sample from distribution.
    std::random_device rd{};
    std::mt19937 gen(rd());
    std::normal_distribution<> norm(mu, sigma);
    double sample = norm(gen);
    return sample;
}
