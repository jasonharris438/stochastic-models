#ifndef STOCHASTIC_MODELS_TRADING_TRADING_LEVELS_EXPONENTIAL_H
#define STOCHASTIC_MODELS_TRADING_TRADING_LEVELS_EXPONENTIAL_H
#include "stochastic_models/trading/exponential_mean_reversion.h"
#include "stochastic_models/trading/optimal_trading.h"
#include "stochastic_models/trading/trading_levels_interface.h"

#include <memory>

/**
 * @brief Class for calculating the optimal trading levels.
 *
 * @param optimizer Pointer to the OptimalTrading instance that is used to
 * find b* - the optimal exit level.
 * @param model Pointer to the stochastic model to use in the function to
 * evluate b*.
 * @param hitting_time_kernel Pointer to the kernel to use in the function
 * to evaluate b*.
 */
class OrnsteinUhlenbeckTradingLevelsExponential : public TradingLevels {
private:
  std::unique_ptr<ExponentialMeanReversion> optimizer;
  std::unique_ptr<StochasticModel> model;
  std::unique_ptr<HittingTimeOrnsteinUhlenbeck> hitting_time_kernel;
  // std::unique_ptr<TradingLevelsParams> params;

public:
  OrnsteinUhlenbeckTradingLevelsExponential(
      const double mu, const double alpha, const double sigma
  );
  const ExponentialMeanReversion* getOptimizer() const;
  const ExponentialMeanReversion* newOptimizer() const;
  const StochasticModel* getModel() const;
  const StochasticModel* newModel() const;
  const HittingTimeOrnsteinUhlenbeck* getHittingTimeKernel() const;
  const HittingTimeOrnsteinUhlenbeck* newHittingTimeKernel() const;
  /**
   * @brief Calculates the hard lower bound constraint for an optimal entry
   * level.
   *
   * This ensures that a reasonable and coherent bound is applied when
   * calculating the optimal entry level.
   *
   * @return const double The hard lower bound constraint for the optimal
   * entry level.
   */
  const double optimalEntryLowerBound() const;
  /**
   * @brief Calculates the hard lower bound constraint for an optimal exit
   * level.
   *
   * The lower bound is the maximum of the cost of trading and the L function.
   * This ensures that a reasonable and coherent bound is applied when
   * calculating the optimal exit level.
   *
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The hard lower bound constraint for the optimal exit
   * level.
   */
  const double optimalExitLowerBound(const double& r, const double& c) const;
  /**
   * @brief Calculates the hard upper bound constraint for an optimal exit
   * level.
   *
   * This ensures that a reasonable and coherent bound is applied when
   * calculating the optimal exit level.
   *
   * @return const double The hard upper bound constraint for the optimal exit
   * level.
   */
  const double optimalExitUpperBound() const;
  /**
   * @brief Calculates the optimal exit level b*
   * for an optimal trading strategy with a stop loss level.
   *
   * @param stop_loss The stop loss level.
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The optimal trading exit level given a stop loss @f[
   b*.
   */
  const double
  optimalExit(const double& stop_loss, const double& r, const double& c) const;
  /**
   * @brief Calculates the optimal exit level b*
   * for an optimal trading strategy.
   *
   * @param optimizer Pointer to the OptimalTrading instance that is used to
   * find b* - the optimal exit level.
   * @param model Pointer to the stochastic model to use in the function to
   * evluate b*.
   * @param hitting_time_kernel Pointer to the kernel to use in the function
   * to evaluate b*.
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The optimal trading exit level b*.
   */
  const double optimalExit(const double& r, const double& c) const;
  /**
   * @brief Calculates the optimal entry level d*
   * for an optimal trading strategy when a stop loss level is provided.
   *
   * @param b_star The optimal exit level b*.
   * @param stop_loss The stop loss level.
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The optimal trading entry level d*.
   */
  const double optimalEntry(
      const double& b_star,
      const double& stop_loss,
      const double& r,
      const double& c
  ) const;
  /**
   * @brief Calculates the lower bound optimal entry level a*
   * for an optimal trading strategy.
   *
   * @param optimizer Pointer to the OptimalTrading instance that is used to
   * find d* - the optimal entry level.
   * @param model Pointer to the stochastic model to use in the function to
   * evluate b*.
   * @param hitting_time_kernel Pointer to the kernel to use in the function
   * to evaluate d*.
   * @param d_star The upper optimal entry level d*.
   * @param b_star The optimal exit level b*.
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The lower bound optimal trading entry level a*
   *.
   */
  const double optimalEntryLower(
      const double& d_star,
      const double& b_star,
      const double& r,
      const double& c
  ) const;
  /**
   * @brief Calculates the lower bound optimal entry level a*
   * for an optimal trading strategy when a stop loss level is provided.
   *
   * @param d_star The upper optimal entry level d*.
   * @param b_star The optimal exit level b*.
   * @param stop_loss The stop loss level.
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The lower bound optimal trading entry level a*
   *.
   */
  const double optimalEntryLower(
      const double& d_star,
      const double& b_star,
      const double& stop_loss,
      const double& r,
      const double& c
  ) const;
  /**
   * @brief Calculates the optimal entry level d*
   * for an optimal trading strategy.
   *
   * @param b_star The optimal exit level b*.
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The optimal trading entry level d*.
   */
  const double
  optimalEntry(const double& b_star, const double& r, const double& c) const;
};

#endif // STOCHASTIC_MODELS_TRADING_TRADING_LEVELS_EXPONENTIAL_H
