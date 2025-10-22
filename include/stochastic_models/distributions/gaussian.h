#ifndef STOCHASTIC_MODELS_DISTRIBUTIONS_GAUSSIAN_H
#define STOCHASTIC_MODELS_DISTRIBUTIONS_GAUSSIAN_H
#include "stochastic_models/distributions/core.h"
/**
 * @file
 * @brief Gaussian (normal) distribution concrete implementation.
 */

/**
 * @brief Gaussian distribution class.
 *
 * Stores mu (mean) and sigma (standard deviation) and implements sampling
 * and CDF helper methods.
 */
class GaussianDistribution : public CoreDistribution {
private:
  // Gaussian distribution mean value.
  const double mu;
  // Gaussian distribution standard deviation value.
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
  /**
   * @brief Helper implementing the error-function-based CDF evaluation.
   * @param x Value at which to evaluate (single-precision is sufficient
   *          for this helper).
   * @return double CDF value.
   */
  double erfGaussianCdf(const float x);

public:
  /**
   * @brief Default constructor producing a standard normal (mu=0, sigma=1).
   */
  GaussianDistribution();
  /**
   * @brief Construct GaussianDistribution with specific mu and sigma.
   * @param mu Mean.
   * @param sigma Standard deviation (must be > 0).
   */
  GaussianDistribution(const double mu, const double sigma);

  /**
   * Returns distribution mean.
   *
   * @returns Const mean of Gaussian Distribution.
   */
  const double getMean() const override;

  /**
   * Returns distribution standard deviation.
   *
   * @returns Const standard deviation of Gaussian Distribution.
   */
  const double getStd() const override;

  /**
   * Uses private erfGaussianCdf method to produce the CDF of the Gaussian
   * distribution function evaluated at the x value provided.
   * Uses mu and sigma parameters set at instantiation.
   *
   * @param x the point at which to evaluate CDF.
   * @returns Gaussian CDF value evaluated at x.
   */
  const double cdf(const double& x) override;

  /**
   * Draws a random sample from normal distribution. Parameterized by
   * mu and sigma private attributes.
   *
   * @param size how many samples to draw. Defaults to 1.
   * @returns Random values drawn from normal distribution.
   */
  std::vector<double> sample(const std::size_t& size = 1) const override;
  ~GaussianDistribution() override;
};
#endif // STOCHASTIC_MODELS_DISTRIBUTIONS_GAUSSIAN_H
