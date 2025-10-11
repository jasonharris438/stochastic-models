#ifndef EXPONENTIAL_MEAN_REVERSION_H
#define EXPONENTIAL_MEAN_REVERSION_H
#include "stochastic_models/sde/stochastic_model.h"
#include "stochastic_models/trading/optimal_trading.h"

/**
 * @brief
 *
 */
class ExponentialMeanReversion : public OptimalTrading {
public:
  /**
   * @brief Construct a new ExponentialMeanReversion object and return on heap
   * memory using the class' copy constructor in the caller instance
   *
   * @return const ExponentialMeanReversion* Pointer to the new instance.
   */
  const ExponentialMeanReversion* clone() const override;
  /**
   * @brief Calculates @f[ F(x;r) @f] in the optimal
   * trading strategy for an exponential mean reverting model.
   *
   * @param hitting_time_kernel Pointer to the kernel instance to use in the
   * function @f[ F(x;r) @f].
   * @param x The current value x to evaluate @f[ F(x;r) @f].
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value of @f[ F(x;r) @f].
   */
  const double
  F(const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates @f[ G(x;r) @f] in the optimal
   * trading strategy for an exponential mean reverting model.
   *
   * @param hitting_time_kernel Pointer to the kernel instance to use in the
   * function @f[ G(x;r) @f].
   * @param x The current value x to evaluate @f[ G(x;r) @f].
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value of @f[ G(x;r) @f].
   */
  const double
  G(const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the value @f[ L^* @f] for the optimal
   * trading strategy for an exponential model.
   *
   * @param hitting_time_kernel Pointer to the kernel instance to
   * use in the function @f[ L^* @f].
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value of @f[ L^* @f].
   */
  const double L_star(
      const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
      const double& r,
      const double& c
  ) const override;
  /**
   * @brief Calculates the value @f[ b @f] for the optimal trading strategy
   * for an an exponential mean reverting model with a stop loss level.
   *
   * @param value The value at which to evaluate the optimal exit level.
   * @param hitting_time_kernel Pointer to the kernel instance to use in the
   * function @f[ b @f].
   * @param stop_loss The stop loss level to use in the optimal trading.
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value @f[ b @f].
   */
  const double
  b(const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& stop_loss,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the value @f[ b @f] in the optimal trading strategy for
   * an exponential mean reverting model.
   *
   * @param value The value at which to evaluate the optimal exit level.
   * @param hitting_time_kernel Pointer to the kernel instance to use in the
   * function @f[ b
   * @f].
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value @f[ b @f].
   */
  const double
  b(const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the value @f[ a @f], the lower bound for the optimal
   * entry region for the optimal trading strategy for an exponential mean
   * reverting model.
   *
   * @param value The value at which to evaluate the optimal entry level.
   * @param hitting_time_kernel Pointer to the kernel instance to use in the
   * function @f[ d
   * @f].
   * @param b_star The optimal exit level calculated by the OptimalTrading
   * instance.
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value @f[ a @f].
   */
  const double
  a(const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star,
    const double& r,
    const double& c) const;
  /**
   * @brief Calculates the value @f[ a @f], the lower bound for the optimal
   * entry region for the optimal trading strategy for an exponential mean
   * reverting model when a stop loss level is provided.
   *
   * @param value The value at which to evaluate the optimal entry level.
   * @param hitting_time_kernel Pointer to the kernel instance to use in the
   * function @f[ d
   * @f].
   * @param b_star The optimal exit level calculated by the OptimalTrading
   * instance.
   * @param stop_loss The stop loss level to use in the optimal trading.
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value @f[ a @f].
   */
  const double
  a(const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star,
    const double& stop_loss,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the value @f[ d @f] for the optimal trading strategy
   * for an exponential mean reverting model when a stop loss level is
   * provided.
   *
   * @param value The value at which to evaluate the optimal entry level.
   * @param hitting_time_kernel Pointer to the kernel instance to use in the
   * function @f[ d
   * @f].
   * @param b_star The optimal exit level calculated by the OptimalTrading
   * instance.
   * @param stop_loss The stop loss level to use in the optimal trading.
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value @f[ d @f].
   */
  const double
  d(const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star,
    const double& stop_loss,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the value @f[ d @f] for the optimal trading strategy
   * for a mean reverting model.
   *
   * @param value The value at which to evaluate the optimal entry level.
   * @param hitting_time_kernel Pointer to the kernel instance to use in the
   * function @f[ d
   * @f].
   * @param b_star The optimal exit level calculated by the OptimalTrading
   * instance.
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value @f[ d @f].
   */
  const double
  d(const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the instantaneous rate of change
   * of the function fn at point x using a hitting time kernel instance.
   *
   * @param fn The function with which to calculate the instantaneous rate of
   * change.
   * @param hitting_time_kernel The mean reverting model instance to use in
   * calculating the instantaneous rate of change.
   * @param x The point at which to calculate the instantaneous rate of
   * change.
   * @param r The discount rate to apply to provide to the function when
   * differentiating.
   * @param c The cost of trading to provide to the function when
   * differentiating.
   * @return const double The instantaneous rate of change of the function at
   * the point x.
   */
  const double instantaneousDifferential(
      ModelFunc fn,
      const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
      const double& x,
      const double& r,
      const double& c
  ) const override;
  /**
   * @brief Calculates the instantaneous rate of change
   * of the function fn at point x using a hitting time kernel instance.
   *
   * @param fn The function with which to calculate the instantaneous rate of
   * change.
   * @param hitting_time_kernel The mean reverting model instance to use in
   * calculating the instantaneous rate of change.
   * @param x The point at which to calculate the instantaneous rate of
   * change.
   * @param b_star The optimal exit level calculated by the OptimalTrading
   * instance.
   * @param r The discount rate to apply to provide to the function when
   * differentiating.
   * @param c The cost of trading to provide to the function when
   * differentiating.
   * @return const double The instantaneous rate of change of the function at
   * the point x.
   */
  const double instantaneousDifferential(
      ModelFunc fn,
      const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
      const double& x,
      const double& b_star,
      const double& r,
      const double& c
  ) const override;
  /**
   * @brief Calculates the instantaneous rate of change
   * of the function fn at point x using a hitting time kernel instance. This
   * overloaded version requires a stop loss level to be provided.
   *
   * @param fn The function with which to calculate the instantaneous rate of
   * change.
   * @param hitting_time_kernel The mean reverting model instance to use in
   * calculating the instantaneous rate of change.
   * @param x The point at which to calculate the instantaneous rate of
   * change.
   * @param b_star The optimal exit level calculated by the OptimalTrading
   * instance.
   * @param stop_loss The stop loss level to use in the optimal trading.
   * @param r The discount rate to apply to provide to the function when
   * differentiating.
   * @param c The cost of trading to provide to the function when
   * differentiating.
   * @return const double The instantaneous rate of change of the function at
   * the point x.
   */
  const double instantaneousDifferential(
      ModelFunc fn,
      const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
      const double& x,
      const double& b_star,
      const double& stop_loss,
      const double& r,
      const double& c
  ) const override;
  /**
   * @brief The value function for the optimal mean reversion trading model.
   *
   * @param hitting_time_kernel The model to use in calculating the value
   * function.
   * @param x The point at which to calculate the instantaneous rate of
   * change.
   * @param b_star The optimal exit level calculated for the model type
   * provided.
   * @param r The discount rate to apply to provide to the function when
   * calculating the value function.
   * @param c The cost of trading to provide to the function when
   * calculating the value function.
   * @return const double The value function evaluated at x given b_star and
   * c.
   */
  const double
  V(const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& b_star,
    const double& r,
    const double& c) const override;
  /**
   * @brief The value function for the optimal trading model for an
   * exponential mean reverting model when a stop loss is provided.
   *
   * @param hitting_time_kernel The model to use in calculating the value
   * function.
   * @param x The point at which to calculate the instantaneous rate of
   * change.
   * @param b_star The optimal exit level calculated for the model type
   * provided.
   * @param stop_loss The stop loss level to use in the optimal trading.
   * @param r The discount rate to apply to provide to the function when
   * calculating the value function.
   * @param c The cost of trading to provide to the function when
   * calculating the value function.
   * @return const double The value function evaluated at x given b_star and
   * c.
   */
  const double
  V(const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& b_star,
    const double& stop_loss,
    const double& r,
    const double& c) const override;
};
#endif // EXPONENTIAL_MEAN_REVERSION_H
