#include <cmath>
#include <random>
#include "gaussian.h"


GaussianDistribution::Gaussian(float mu, float sigma) :
    mu(mu), sigma(sigma) { }

GaussianDistribution::Gaussian() :
    GaussianDistribution(0, 1.0) { }

float GaussianDistribution::erfGaussianCdf(float x){    // Core calculation for cdf evaluated at x.
    return (1.0/2.0) * (1.0 + std::erf(x));
}

float GaussianDistribution::Cdf(float x){    // Produces cdf for a given x.
    float val = (x - mu) / (sigma * std::sqrt(2));
    return GaussianDistribution::erfGaussianCdf(val);
}

float GaussianDistribution::Sample(std::mt19937 &gen){   // Draws random sample from distribution.
    std::normal_distribution<float> norm(mu, sigma);
    return norm;
}
