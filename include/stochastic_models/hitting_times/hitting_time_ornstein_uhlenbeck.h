#ifndef STOCHASTIC_MODELS_HITTING_TIMES_HITTING_TIME_ORNSTEIN_UHLENBECK_H
#define STOCHASTIC_MODELS_HITTING_TIMES_HITTING_TIME_ORNSTEIN_UHLENBECK_H
#include "stochastic_models/sde/ornstein_uhlenbeck.h"

/**
 * @file
 * @brief Kernel used when computing OU hitting-time densities.
 */

/**
 * @brief Hitting time kernel that implements the integrand and helper
 * functions specific to the O-U model that are used by the hitting-time numeric
 * routines.
 */
class HittingTimeOrnsteinUhlenbeck {
private:
  const double mu;
  const double alpha;
  const double sigma;

public:
  HittingTimeOrnsteinUhlenbeck(
      const double mu, const double alpha, const double sigma
  );
  HittingTimeOrnsteinUhlenbeck(const HittingTimeOrnsteinUhlenbeck& other);
  /**
   * @brief Return a heap-allocated copy of this instance.
   *
   * @return const HittingTimeOrnsteinUhlenbeck* Pointer to the new instance.
   */
  const HittingTimeOrnsteinUhlenbeck* clone() const;
  /**
   * @brief Evaluate the integrand core S(x) of the hitting-time density.
   *
   * @param x The point at which to evaluate the hitting time density function
   * core.
   * @return const double The hitting time density function core evaluated at
   * x.
   */
  const double hittingTimeDensityCore(const double& x) const;
  /**
   * @brief Kernel function F(x,u,r) used in optimal trading integrals.
   * @param x The point x at which to evaluate the first hitting time density
   * function kernel.
   * @param u The point u at which to evaluate the first hitting time density
   * function kernel.
   * @param r The value r indicating the discount rate.
   * @return const double The first hitting time density function kernel
   * evaluated at x.
   */
  const double
  optimalTradingFCore(const double& x, const double& u, const double& r) const;
  /**
   * @brief Kernel function G(x,u,r) used in optimal trading integrals.
   * @param x The point x at which to evaluate the first hitting time density
   * function kernel.
   * @param u The point u at which to evaluate the first hitting time density
   * function kernel.
   * @param r The value r indicating the discount rate.
   * @return const double The first hitting time density function kernel
   * evaluated at x.
   */
  const double
  optimalTradingGCore(const double& x, const double& u, const double& r) const;
  /**
   * @brief Computes the L*(r,c) optimal trading helper function.
   *
   * @param r The value r indicating the discount rate.
   * @param c The value c indicating the cost of trading.
   * @return const double The optimal trading funcion L*.
   */
  const double optimalTradingLCore(const double& r, const double& c) const;
  ~HittingTimeOrnsteinUhlenbeck() = default;
};
#endif // STOCHASTIC_MODELS_HITTING_TIMES_HITTING_TIME_ORNSTEIN_UHLENBECK_H
