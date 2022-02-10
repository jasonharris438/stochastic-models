#include <cmath>
#include <random>
#include "gaussian.h"


GaussianDistribution::GaussianDistribution(const double mu, const double sigma) :
    mu(mu), sigma(sigma) { }

GaussianDistribution::GaussianDistribution() :
    GaussianDistribution(0, 1.0) { }

const double &const GaussianDistribution::getMean(){    // Returns const reference to class mean.
    return mu;
}

const double &const GaussianDistribution::getStd(){    // Returns const reference to class std.
    return sigma;
}

double GaussianDistribution::erfGaussianCdf(const float x){    // Core calculation for cdf evaluated at x.
    return (1.0/2.0) * (1.0 + std::erf(x));
}

double GaussianDistribution::Cdf(const double x){    // Produces cdf for a given x.
    const double val = (x - mu) / (sigma * std::sqrt(2));
    return GaussianDistribution::erfGaussianCdf(val);
}

double GaussianDistribution::Sample(std::mt19937 &gen) const {   // Draws random sample from distribution.
    std::normal_distribution<> norm(mu, sigma);
    return norm(gen);
}
