#ifndef TRADING_LEVELS_H
#define TRADING_LEVELS_H
#include "stochastic_models/trading/optimal_mean_reversion.h"
#include "stochastic_models/trading/optimal_trading.h"
#include "stochastic_models/trading/trading_levels_interface.h"

#include <memory>

/**
 * @brief Class for calculating the optimal trading levels.
 *
 * @param optimizer Pointer to the OptimalTrading instance that is used to
 * find @f[ b* @f] - the optimal exit level.
 * @param model Pointer to the stochastic model to use in the function to
 * evluate @f[ b* @f].
 * @param hitting_time_kernel Pointer to the kernel to use in the function
 * to evaluate @f[ b* @f].
 */
class OrnsteinUhlenbeckTradingLevels : public TradingLevels {
private:
  std::unique_ptr<OptimalMeanReversion> optimizer;
  std::unique_ptr<StochasticModel> model;
  std::unique_ptr<HittingTimeOrnsteinUhlenbeck> hitting_time_kernel;
  // std::unique_ptr<TradingLevelsParams> params;

public:
  OrnsteinUhlenbeckTradingLevels(
      const double mu, const double alpha, const double sigma
  );
  const OptimalMeanReversion* getOptimizer() const;
  const OptimalMeanReversion* newOptimizer() const;
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
   * @brief Calculates the optimal exit level @f[ b* @f]
   * for an optimal trading strategy with a stop loss level.
   *
   * @param stop_loss The stop loss level.
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The optimal trading exit level given a stop loss @f[
   b* @f].
   */
  const double
  optimalExit(const double& stop_loss, const double& r, const double& c) const;
  /**
   * @brief Calculates the optimal exit level @f[ b* @f]
   * for an optimal trading strategy.
   *
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The optimal trading exit level @f[ b* @f].
   */
  const double optimalExit(const double& r, const double& c) const;
  /**
   * @brief Calculates the optimal entry level @f[ d* @f]
   * for an optimal trading strategy when a stop loss level is provided.
   *
   * @param b_star The optimal exit level @f[ b* @f].
   * @param stop_loss The stop loss level.
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The optimal trading entry level @f[ d* @f].
   */
  const double optimalEntry(
      const double& b_star,
      const double& stop_loss,
      const double& r,
      const double& c
  ) const;
  /**
   * @brief Calculates the lower bound optimal entry level @f[ a* @f]
   * for an optimal trading strategy.
   *
   * @param d_star The upper optimal entry level @f[ d* @f].
   * @param b_star The optimal exit level @f[ b* @f].
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The lower bound optimal trading entry level @f[ a*
   * @f].
   */
  const double optimalEntryLower(
      const double& d_star,
      const double& b_star,
      const double& r,
      const double& c
  ) const;
  /**
   * @brief Calculates the lower bound optimal entry level @f[ a* @f]
   * for an optimal trading strategy when a stop loss level is provided.
   *
   * @param d_star The upper optimal entry level @f[ d* @f].
   * @param b_star The optimal exit level @f[ b* @f].
   * @param stop_loss The stop loss level.
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The lower bound optimal trading entry level @f[ a*
   * @f].
   */
  const double optimalEntryLower(
      const double& d_star,
      const double& b_star,
      const double& stop_loss,
      const double& r,
      const double& c
  ) const;
  /**
   * @brief Calculates the optimal entry level @f[ d* @f]
   * for an optimal trading strategy.
   *
   * @param b_star The optimal exit level @f[ b* @f].
   * @param r The discount rate to apply to the optimal mean reversion trading
   * problem.
   * @param c The cost of trading.
   * @return const double The optimal trading entry level @f[ d* @f].
   */
  const double
  optimalEntry(const double& b_star, const double& r, const double& c) const;
};

#endif // TRADING_LEVELS_H
