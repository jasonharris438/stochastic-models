#ifndef STOCHASTIC_MODELS_DISTRIBUTIONS_CORE_H
#define STOCHASTIC_MODELS_DISTRIBUTIONS_CORE_H
#include <vector>

/**
 * @file
 * @brief Core distribution abstractions used by concrete distribution classes.
 */

/**
 * @brief Abstract base class for statistical distributions used in the library.
 *
 * Concrete distributions must derive from this interface and provide mean/std,
 * sampling and CDF computation.
 */
class CoreDistribution {
public:
  /**
   * @brief Virtual destructor for interface cleanup.
   */
  virtual ~CoreDistribution() = 0;
  /**
   * @brief Return the distribution mean.
   * @return const double Mean value.
   */
  virtual const double getMean() const = 0;

  /**
   * @brief Return the distribution standard deviation.
   * @return const double Standard deviation.
   */
  virtual const double getStd() const = 0;

  /**
   * @brief Compute the cumulative distribution function at x.
   *
   * @param x Point at which to evaluate the CDF.
   * @return const double CDF(x).
   */
  virtual const double cdf(const double& x) = 0;

  /**
   * @brief Draw random samples from the distribution.
   *
   * @param size Number of samples to draw (default 1).
   * @return std::vector<double> Vector of drawn samples.
   */
  virtual std::vector<double> sample(const std::size_t& size = 1) const = 0;
};
#endif // STOCHASTIC_MODELS_DISTRIBUTIONS_CORE_H
