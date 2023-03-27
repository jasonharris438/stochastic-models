#ifndef _CORE_H_
#define _CORE_H_
#include <vector>
/**
 * Base core distribution class.
 */
class CoreDistribution {
   public:
    /**
     * Returns const reference to distribution mean.
     *
     * @returns Const reference to mean of Gaussian Distribution.
     */
    virtual const double &getMean() const = 0;

    /**
     * Returns const reference to distribution standard deviation.
     *
     * @returns Const reference to standard deviation of Gaussian Distribution.
     */
    virtual const double &getStd() const = 0;

    /**
     * Uses private erfGaussianCdf method to produce the CDF of the Gaussian
     * distribution function evaluated at the x value provided.
     * Uses mu and sigma parameters set at instantiation.
     *
     * @param x the point at which to evaluate CDF.
     * @returns Gaussian CDF value evaluated at x.
     */
    virtual const double Cdf(const double &x) = 0;

    /**
     * Draws a random sample from normal distribution. Parameterized by
     * mu and sigma private attributes.
     *
     * @param size how many samples to draw. Defaults to 1.
     * @returns Random values drawn from normal distribution.
     */
    virtual std::vector<double> Sample(std::size_t size = 1) const = 0;
    virtual ~CoreDistribution() = 0;
};
#endif  // _GAUSSIAN_H_
