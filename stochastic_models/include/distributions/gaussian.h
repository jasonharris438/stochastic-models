#ifndef _GAUSSIAN_H_
#define _GAUSSIAN_H_
#include "core.h"

/**
 * Gaussian distribution class. Takes two floats as inputs "mu" and "sigma" for
 * mean and standard deviation and implements several functions related
 * to probabilistic characteristics of the Gaussian distribution.
 */
class GaussianDistribution : public CoreDistribution {
   private:
    const double mu;
    const double sigma;

    /**
     * Uses Error Function to produce the CDF of the Gaussian
     * distribution function evaluated at the x value provided.
     * Uses mu and sigma parameters set at instantiation.
     *
     * article:
     * https://stats.stackexchange.com/questions/187828/how-are-the-error-function-and-standard-normal-distribution-function-related
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
    const double &getMean() const override;

    /**
     * Returns const reference to distribution standard deviation.
     *
     * @returns Const reference to standard deviation of Gaussian Distribution.
     */
    const double &getStd() const override;

    /**
     * Uses private erfGaussianCdf method to produce the CDF of the Gaussian
     * distribution function evaluated at the x value provided.
     * Uses mu and sigma parameters set at instantiation.
     *
     * @param x the point at which to evaluate CDF.
     * @returns Gaussian CDF value evaluated at x.
     */
    double Cdf(const double x) override;

    /**
     * Draws a random sample from normal distribution. Parameterized by
     * mu and sigma private attributes.
     *
     * @param size how many samples to draw. Defaults to 1.
     * @returns Random values drawn from normal distribution.
     */
    std::vector<double> Sample(std::size_t size = 1) const override;
    ~GaussianDistribution(){};
};

#endif  // _GAUSSIAN_H_
