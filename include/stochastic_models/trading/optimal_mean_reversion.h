#ifndef STOCHASTIC_MODELS_TRADING_OPTIMAL_MEAN_REVERSION_H
#define STOCHASTIC_MODELS_TRADING_OPTIMAL_MEAN_REVERSION_H
#include "stochastic_models/hitting_times/hitting_time_ornstein_uhlenbeck.h"
#include "stochastic_models/sde/stochastic_model.h"
#include "stochastic_models/trading/optimal_trading.h"

/**
 * @brief Optimal mean reversion trading model parameters.
 * @param hitting_time_kernel The hitting time kernel instance to use in
 * numerical integration.
 * @param x The current value of the stochastic model.
 * @param r The discount rate to apply to the optimal trading problem.
 *
 */
struct OptimalMeanReversionParams {
  const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel;
  const double& x;
  const double& r;

  ~OptimalMeanReversionParams();
};
/**
 * @brief Used when evaluating the optimal mean reversion trading
 * value function V(x).
 *
 * @param x The point at which to evaluate the function V(x).
 * @param params The parameters provided to use to implement the function @f[
 * V(x) with respect to the model type.
 * @return double The value of the integral V(x) at the point x.
 */
double valueFunc(double x, void* params);
/**
 * @brief Used when evaluating the optimal mean reversion trading
 * value function V(x) when a stop loss is provided.
 *
 * @param x The point at which to evaluate the function V(x).
 * @param params The parameters provided to use to implement the function @f[
 * V(x) with respect to the model type.
 * @return double The value of the integral V(x) at the point x.
 */
double valueFuncStopLoss(double x, void* params);
/**
 * @brief Used when integrating to evaluate the optimal mean reversion trading
 * model function G(x;r).
 *
 * @param x The point at which to evaluate the function G(x;r).
 * @param params The parameters provided to use to implement the function @f[
 * G(x;r) with respect to the model type.
 * @return double The value of the integral G(x;r) at the point x.
 */
double funcIntegrateG(double x, void* params);
/**
 * @brief Used when integrating to evaluate the optimal mean reversion trading
 * model function F(x;r).
 *
 * @param x The point at which to evaluate the function F(x;r).
 * @param params The parameters provided to use to implement the function @f[
 * F(x;r) with respect to the model type.
 * @return double The value of the integral F(x;r) at the point x.
 */
double funcIntegrateF(double x, void* params);
/**
 * @brief Used when calculating the optimal mean reversion trading model
 * function F(x;r).
 *
 * @param x The point at which to evaluate the antiderivative.
 * @param params The parameters provided via the numerical integration function
 * that contains the model and discount rate to use.
 * @return double The value of the integral F(x;r) at the point x.
 */
double funcOptimalMeanReversionF(double x, void* params);
/**
 * @brief Used when calculating the optimal mean reversion trading model
 * function G(x;r).
 *
 * @param x The point at which to evaluate the antiderivative.
 * @param params The parameters provided via the numerical integration function
 * that contains the model and discount rate to use.
 * @return double The value of the integral G(x;r) at the point x.
 */
double funcOptimalMeanReversionG(double x, void* params);

/**
 * @brief Concrete class that implements the optimal trading strategy for a mean
 * reverting model.
 */
class OptimalMeanReversion : public OptimalTrading {
public:
  /**
   * @brief Construct a new OptimalMeanReversion object and return on heap
   * memory using the class' copy constructor in the caller instance
   *
   * @return const OptimalMeanReversion* Pointer to the new instance.
   */
  const OptimalMeanReversion* clone() const override;
  /**
   * @brief Calculates F(x;r) in the optimal
   * trading strategy for a mean reverting model.
   *
   * @param hitting_time_kernel Pointer to the hitting time kernel instance to
   * use in the function F(x;r).
   * @param x The current value x to evaluate F(x;r).
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value of F(x;r).
   */
  const double
  F(const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates G(x;r) in the optimal
   * trading strategy for a mean reverting model.
   *
   * @param hitting_time_kernel Pointer to the hitting time kernel instance to
   * use in the function G(x;r).
   * @param x The current value x to evaluate G(x;r).
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value of G(x;r).
   */
  const double
  G(const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the value L* for the optimal trading strategy.
   *
   * @param hitting_time_kernel Pointer to the hitting time kernel to use in
   * the function
   * L*.
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value of L*.
   */
  const double L_star(
      const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
      const double& r,
      const double& c
  ) const override;
  /**
   * @brief Calculates the value b for the optimal trading strategy
   * for a mean reverting model with a stop loss level.
   *
   * @param value The value at which to evaluate the optimal exit level.
   * @param hitting_time_kernel Pointer to the hitting time kernel to use in
   * the function b
   *.
   * @param stop_loss The stop loss level to use in the optimal trading.
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value b.
   */
  const double
  b(const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& stop_loss,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the value b in the optimal trading strategy for
   * a mean reverting model.
   *
   * @param value The value at which to evaluate the optimal exit level.
   * @param hitting_time_kernel Pointer to the hitting time kernel to use in
   * the function b
   *.
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value b.
   */
  const double
  b(const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the value a, the lower bound for the optimal
   * entry region for the optimal trading strategy for a mean reverting model.
   *
   * @param value The value at which to evaluate the optimal entry level.
   * @param hitting_time_kernel Pointer to the hitting time kernel to use in
   * the function d
   *.
   * @param b_star The optimal exit level calculated by the OptimalTrading
   * instance.
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value a.
   */
  const double
  a(const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the value a, the lower bound for the optimal
   * entry region for the optimal trading strategy for a mean reverting model
   * when a stop loss level is provided.
   *
   * @param value The value at which to evaluate the optimal entry level.
   * @param hitting_time_kernel Pointer to the hitting time kernel to use in
   * the function d
   *.
   * @param b_star The optimal exit level calculated by the OptimalTrading
   * instance.
   * @param stop_loss The stop loss level to use in the optimal trading.
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value a.
   */
  const double
  a(const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star,
    const double& stop_loss,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the value d for the optimal trading strategy
   * for a mean reverting model when a stop loss level is provided.
   *
   * @param value The value at which to evaluate the optimal entry level.
   * @param hitting_time_kernel Pointer to the hitting time kernel to use in
   * the function d
   *.
   * @param b_star The optimal exit level calculated by the OptimalTrading
   * instance.
   * @param stop_loss The stop loss level to use in the optimal trading.
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value d.
   */
  const double
  d(const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star,
    const double& stop_loss,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the value d for the optimal trading strategy
   * for a mean reverting model.
   *
   * @param value The value at which to evaluate the optimal entry level.
   * @param hitting_time_kernel Pointer to the hitting time kernel to use in
   * the function d
   *.
   * @param b_star The optimal exit level calculated by the OptimalTrading
   * instance.
   * @param r The discount rate to apply to the optimal trading problem.
   * @param c The cost of trading.
   * @return const double The value d.
   */
  const double
  d(const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star,
    const double& r,
    const double& c) const override;
  /**
   * @brief Calculates the instantaneous rate of change
   * of the function fn at point x using a mean reverting model instance
   * pointed to by model.
   *
   * @param fn The function with which to calculate the instantaneous rate of
   * change.
   * @param hitting_time_kernel The kernel instance to use in calculating the
   * instantaneous rate of change.
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
   * of the function fn at point x using a hitting_time_kernel instance.
   *
   * @param fn The function with which to calculate the instantaneous rate of
   * change.
   * @param hitting_time_kernel The kernel instance to use in calculating the
   * instantaneous rate of change.
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
   * of the function fn at point x using a hitting_time_kernel instance. This
   * overloaded version requires a stop loss level to be provided.
   *
   * @param fn The function with which to calculate the instantaneous rate of
   * change.
   * @param hitting_time_kernel The kernel instance to use in calculating the
   * instantaneous rate of change.
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
   * @param hitting_time_kernel The kernel instance to use in calculating the
   * value function.
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
   * @brief The value function for the optimal trading model when a stop loss
   * is provided.
   *
   * @param hitting_time_kernel The kernel instance to use in calculating the
   * value function.
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
#endif // STOCHASTIC_MODELS_TRADING_OPTIMAL_MEAN_REVERSION_H
