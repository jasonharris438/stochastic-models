#ifndef STOCHASTIC_MODELS_TRADING_TRADING_LEVELS_H
#define STOCHASTIC_MODELS_TRADING_TRADING_LEVELS_H
#include "stochastic_models/trading/optimal_mean_reversion.h"
#include "stochastic_models/trading/optimal_trading.h"
#include "stochastic_models/trading/trading_levels_interface.h"

#include <memory>

/**
 * @file
 * @brief High-level helpers to compute optimal trading entry/exit levels.
 */

/**
 * @brief Class for calculating optimal trading levels using the
 * Ornstein-Uhlenbeck process.
 *
 * This class implements trading level calculations based on the
 * Ornstein-Uhlenbeck stochastic process. It provides methods to calculate
 * optimal entry and exit levels for trading strategies with optional stop-loss
 * mechanisms.
 *
 * Based on Leung, T., & Li, X. (2015). Optimal mean reversion trading book.
 */
class OrnsteinUhlenbeckTradingLevels : public TradingLevels {
private:
  /**
   * @brief Unique pointer to the optimal mean reversion optimizer used for
   *              calculating trading strategy parameters
   */
  std::unique_ptr<OptimalMeanReversion> optimizer;
  /**
   * @brief Unique pointer to the stochastic model implementation representing
   *              the Ornstein-Uhlenbeck process
   */
  std::unique_ptr<StochasticModel> model;
  /**
   * @brief Unique pointer to the hitting time kernel specifically
   *        designed for Ornstein-Uhlenbeck processes, used for
   *        calculating expected hitting times between levels
   */
  std::unique_ptr<HittingTimeOrnsteinUhlenbeck> hitting_time_kernel;

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

#endif // STOCHASTIC_MODELS_TRADING_TRADING_LEVELS_H
