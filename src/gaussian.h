#ifndef Gaussian_H_
#define Gaussian_H_

#include <random>

/**
 * Gaussian distribution class. Takes two floats as inputs "mu" and "sigma" for 
 * mean and standard deviation and implements several functions related
 * to probabilistic characteristics of the Gaussian distribution.
 */
class GaussianDistribution {
private:
    float mu;
    float sigma;

    /**
     * Uses Error Function to produce the CDF of the Gaussian
     * distribution function evaluated at the x value provided. 
     * Uses mu and sigma parameters set at instantiation.
     * 
     * article: https://stats.stackexchange.com/questions/187828/how-are-the-error-function-and-standard-normal-distribution-function-related
     * @param x the point at which to evaluate CDF.
     * @returns Gaussian CDF value evaluated at x.
     */
    float erfGaussianCdf(float x);

public:
    GaussianDistribution();
    GaussianDistribution(float mu, float sigma);

    /**
     * Uses private erfGaussianCdf method to produce the CDF of the Gaussian
     * distribution function evaluated at the x value provided. 
     * Uses mu and sigma parameters set at instantiation.
     * 
     * @param x the point at which to evaluate CDF.
     * @returns Gaussian CDF value evaluated at x.
     */
    float Cdf(float x);

    /**
     * Draws a random sample from normal distribution. Parameterized by
     * mu and sigma private attributes.
     * 
     * @param gen Reference to mersenne_twister_engine object passed to be
     * used in random number generation.
     * @returns Single random value drawn from normal distribution.
     */
    float Sample(std::mt19937 &gen);
};

#endif