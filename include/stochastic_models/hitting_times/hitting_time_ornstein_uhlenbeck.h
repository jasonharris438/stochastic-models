#ifndef HITTING_TIME_ORNSTEIN_UHLENBECK_H
#define HITTING_TIME_ORNSTEIN_UHLENBECK_H
#include "stochastic_models/sde/ornstein_uhlenbeck.h"
/**
 * Hitting time kernel class that handles internals for calculating the
 * hitting time density of trading levels.
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
   * @brief Construct a new HittingTimeOrnsteinUhlenbeck object and return on
   * heap memory using the class' copy constructor in the caller instance
   *
   * @return const HittingTimeOrnsteinUhlenbeck* Pointer to the new instance.
   */
  const HittingTimeOrnsteinUhlenbeck* clone() const;
  /**
   * @brief Implements the core components of the hitting time density
   * function for the model.
   *
   * @param x The point at which to evaluate the hitting time density function
   * core.
   * @return const double The hitting time density function core evaluated at
   * x.
   */
  const double hittingTimeDensityCore(const double& x) const;
  /**
   * @brief Implements the kernel of the first hitting time density function
   * @f[ F(x) @f].
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
   * @brief Implements the kernel of the first hitting time density function
   * @f[ G(x) @f].
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
   * @brief Implements the kernel of the optimal trading function @f[ L^* @f].
   *
   * @param r The value r indicating the discount rate.
   * @param c The value c indicating the cost of trading.
   * @return const double The optimal trading funcion @f[ L^* @f].
   */
  const double optimalTradingLCore(const double& r, const double& c) const;
  ~HittingTimeOrnsteinUhlenbeck() = default;
};
#endif // HITTING_TIME_ORNSTEIN_UHLENBECK_H
