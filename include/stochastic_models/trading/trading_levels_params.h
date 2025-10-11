#ifndef TRADING_LEVELS_PARAMS_H
#define TRADING_LEVELS_PARAMS_H
#include "stochastic_models/trading/optimal_trading.h"

/**
 * @brief Optimal mean reversion trading model parameters for finding the
 * optimal exit level when a stop loss is provided.
 * @param optimizer The optimal trading instance to use when calculating the
 * optimal exit level.
 * @param hitting_time_kernel The kernel to use when calculating the optimal
 * exit level.
 * @param stop_loss The stop loss level.
 * @param r The discount rate to apply to the optimal trading problem.
 * @param c The cost of trading.
 */
struct ExitLevelStopLossParams {
  const OptimalTrading* optimizer;
  const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel;
  const double& stop_loss;
  const double& r;
  const double& c;

  ~ExitLevelStopLossParams();
};
/**
 * @brief Optimal mean reversion trading model parameters for finding the
 * optimal exit level.
 *
 * @param optimizer The optimal trading instance to use when calculating the
 * optimal exit level.
 * @param hitting_time_kernel The kernel to use when calculating the optimal
 * exit level.
 * @param r The discount rate to apply to the optimal trading problem.
 * @param c The cost of trading.
 */
struct ExitLevelParams {
  const OptimalTrading* optimizer;
  const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel;
  const double& r;
  const double& c;

  ~ExitLevelParams();
};
/**
 * @brief Optimal mean reversion trading model parameters for finding the
 * optimal entry level when a stop loss is provided.
 *
 * @param optimizer The optimal trading instance to use when calculating optimal
 * optimal entry level.
 * @param hitting_time_kernel The kernel to use when calculating the optimal
 * entry level.
 * @param b_star The optimal exit level.
 * @param stop_loss The stop loss level.
 * @param r The discount rate to apply to the optimal trading problem.
 * @param c The cost of trading.
 */
struct EntryLevelStopLossParams {
  const OptimalTrading* optimizer;
  const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel;
  const double& b_star;
  const double& stop_loss;
  const double& r;
  const double& c;

  ~EntryLevelStopLossParams();
};
/**
 * @brief Optimal mean reversion trading model parameters for finding the
 * optimal entry level.
 *
 * @param optimizer The optimal trading instance to use when calculating optimal
 * optimal entry level.
 * @param hitting_time_kernel The kernel to use when calculating the optimal
 * entry level.
 * @param b_star The optimal exit level.
 * @param r The discount rate to apply to the optimal trading problem.
 * @param c The cost of trading.
 */
struct EntryLevelParams {
  const OptimalTrading* optimizer;
  const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel;
  const double& b_star;
  const double& r;
  const double& c;

  ~EntryLevelParams();
};
/**
 * @brief Function to evaluate the optimal entry level function @f[ a @f] in the
 * optimal trading problem.
 *
 * @param x The value at which to evaluate the optimal entry level function.
 * @param params Pointer to the EntryLevelParams struct containing the
 * parameters used in evaluating the optimal entry level function.
 * @return double The value of the optimal entry level function @f[ d @f] at the
 * point x.
 */
double funcOptimalMeanReversionA(double x, void* params);
/**
 * @brief Function to evaluate the optimal entry level function @f[ a @f] in the
 * optimal trading problem when a stop loss is provided.
 *
 * @param x The value at which to evaluate the optimal entry level function.
 * @param params Pointer to the EntryLevelStopLossParams struct containing the
 * parameters used in evaluating the optimal entry level function.
 * @return double The value of the optimal entry level function @f[ d @f] at the
 * point x.
 */
double funcOptimalMeanReversionStopLossA(double x, void* params);
/**
 * @brief Function to evaluate the optimal entry level function @f[ d @f] in the
 * optimal trading problem when a stop loss is provided.
 *
 * @param x The value at which to evaluate the optimal entry level function.
 * @param params Pointer to the EntryLevelStopLossParams struct containing the
 * parameters used in evaluating the optimal entry level function.
 * @return double The value of the optimal entry level function @f[ d @f] at the
 * point x.
 */
double funcOptimalMeanReversionStopLossD(double x, void* params);
/**
 * @brief Function to evaluate the optimal entry level function @f[ d @f] in the
 * optimal trading problem.
 *
 * @param x The value at which to evaluate the optimal entry level function.
 * @param params Pointer to the EntryLevelParams struct containing the
 * parameters used in evaluating the optimal entry level function.
 * @return double The value of the optimal entry level function @f[ d @f] at the
 * point x.
 */
double funcOptimalMeanReversionD(double x, void* params);
/**
 * @brief Function to evaluate the optimal exit level function @f[ b @f] in the
 * optimal trading problem when there is a stop loss level provided.
 *
 * @param x The value at which to evaluate the optimal exit level function.
 * @param params Pointer to the ExitLevelStopLossParams struct containing the
 * parameters used in evaluating the optimal exit level function.
 * @return double The value of the optimal exit level function @f[ b @f] at the
 * point x.
 */
double funcOptimalMeanReversionStopLossB(double x, void* params);
/**
 * @brief Function to evaluate the optimal exit level function @f[ b @f] in the
 * optimal trading problem when a stop loss is provided.
 *
 * @param x The value at which to evaluate the optimal exit level function.
 * @param params Pointer to the ExitLevelParams struct containing the
 * parameters used in evaluating the optimal exit level function.
 * @return double The value of the optimal exit level function @f[ b @f] at the
 * point x.
 */
double funcOptimalMeanReversionB(double x, void* params);
#endif // TRADING_LEVELS_PARAMS_H
