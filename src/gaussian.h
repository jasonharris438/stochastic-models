#ifndef _GAUSSIAN_H_
#define _GAUSSIAN_H_
#include <random>

/**
 * Gaussian distribution class. Takes two floats as inputs "mu" and "sigma" for 
 * mean and standard deviation and implements several functions related
 * to probabilistic characteristics of the Gaussian distribution.
 */
class GaussianDistribution {
private:
    const double mu;
    const double sigma;

    /**
     * Uses Error Function to produce the CDF of the Gaussian
     * distribution function evaluated at the x value provided. 
     * Uses mu and sigma parameters set at instantiation.
     * 
     * article: https://stats.stackexchange.com/questions/187828/how-are-the-error-function-and-standard-normal-distribution-function-related
     * @param x the point at which to evaluate CDF.
     * @returns Gaussian CDF value evaluated at x.
     */
    double erfGaussianCdf(const float x);

public:
    GaussianDistribution();
    GaussianDistribution(const double mu, const double sigma);

    /**
     * Returns const reference to distribution mean.
     * 
     * @returns Const reference to mean of Gaussian Distribution.
     */
    const double &getMean();

    /**
     * Returns const reference to distribution standard deviation.
     * 
     * @returns Const reference to standard deviation of Gaussian Distribution.
     */
    const double &getStd();

    /**
     * Uses private erfGaussianCdf method to produce the CDF of the Gaussian
     * distribution function evaluated at the x value provided. 
     * Uses mu and sigma parameters set at instantiation.
     * 
     * @param x the point at which to evaluate CDF.
     * @returns Gaussian CDF value evaluated at x.
     */
    double Cdf(const double x);

    /**
     * Draws a random sample from normal distribution. Parameterized by
     * mu and sigma private attributes.
     * 
     * @param gen Reference to mersenne_twister_engine object passed to be
     * used in random number generation.
     * @returns Single random value drawn from normal distribution.
     */
    double Sample(std::mt19937 &gen) const;
};

#endif // _GAUSSIAN_H_
