#include <cmath>
#include <random>
#include "gaussian.h"


GaussianDistribution::GaussianDistribution(const float mu, const float sigma) :
    mu(mu), sigma(sigma) { }

GaussianDistribution::GaussianDistribution() :
    GaussianDistribution(0, 1.0) { }

const float *GaussianDistribution::getMean(){    // Returns pointer to class mean.
    const float *ptr = &mu;
    return ptr;
}

const float *GaussianDistribution::getStd(){    // Returns pointer to class std.
    const float *ptr = &sigma;
    return ptr;
}

float GaussianDistribution::erfGaussianCdf(const float x){    // Core calculation for cdf evaluated at x.
    return (1.0/2.0) * (1.0 + std::erf(x));
}

float GaussianDistribution::Cdf(float x){    // Produces cdf for a given x.
    const float val = (x - mu) / (sigma * std::sqrt(2));
    return GaussianDistribution::erfGaussianCdf(val);
}

double GaussianDistribution::Sample(std::mt19937 &gen) const {   // Draws random sample from distribution.
    std::normal_distribution<> norm(mu, sigma);
    return norm(gen);
}
