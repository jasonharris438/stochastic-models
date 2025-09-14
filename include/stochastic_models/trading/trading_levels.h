#ifndef TRADING_LEVELS_H
#define TRADING_LEVELS_H
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
/**
 * @brief Class for calculating the optimal trading levels.
 *
 */
class TradingLevels {
   public:
    /**
     * @brief Calculates the optimal exit level @f[ b* @f]
     * for an optimal trading strategy with a stop loss level.
     *
     * @param optimizer Pointer to the OptimalTrading instance that is used to
     * find @f[ b* @f] - the optimal exit level.
     * @param model Pointer to the stochastic model to use in the function to
     * evluate @f[ b* @f].
     * @param hitting_time_kernel Pointer to the kernel to use in the function
     * to evaluate @f[ b* @f].
     * @param stop_loss The stop loss level.
     * @param r The discount rate to apply to the optimal mean reversion trading
     * problem.
     * @param c The cost of trading.
     * @return const double The optimal trading exit level given a stop loss @f[
     b* @f].
     */
    const double optimalExit(
        const OptimalTrading* optimizer, const StochasticModel* model,
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& stop_loss, const double& r, const double& c) const;
    /**
     * @brief Calculates the optimal exit level @f[ b* @f]
     * for an optimal trading strategy.
     *
     * @param optimizer Pointer to the OptimalTrading instance that is used to
     * find @f[ b* @f] - the optimal exit level.
     * @param model Pointer to the stochastic model to use in the function to
     * evluate @f[ b* @f].
     * @param hitting_time_kernel Pointer to the kernel to use in the function
     * to evaluate @f[ b* @f].
     * @param r The discount rate to apply to the optimal mean reversion trading
     * problem.
     * @param c The cost of trading.
     * @return const double The optimal trading exit level @f[ b* @f].
     */
    const double optimalExit(
        const OptimalTrading* optimizer, const StochasticModel* model,
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& r, const double& c) const;
    /**
     * @brief Calculates the optimal entry level @f[ d* @f]
     * for an optimal trading strategy when a stop loss level is provided.
     *
     * @param optimizer Pointer to the OptimalTrading instance that is used to
     * find @f[ d* @f] - the optimal entry level.
     * @param model Pointer to the stochastic model to use in the function to
     * evluate @f[ b* @f].
     * @param hitting_time_kernel Pointer to the kernel to use in the function
     * to evaluate @f[ d* @f].
     * @param b_star The optimal exit level @f[ b* @f].
     * @param stop_loss The stop loss level.
     * @param r The discount rate to apply to the optimal mean reversion trading
     * problem.
     * @param c The cost of trading.
     * @return const double The optimal trading entry level @f[ d* @f].
     */
    const double optimalEntry(
        const OptimalTrading* optimizer, const StochasticModel* model,
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& b_star, const double& stop_loss, const double& r,
        const double& c) const;
    /**
     * @brief Calculates the lower bound optimal entry level @f[ a* @f]
     * for an optimal trading strategy.
     *
     * @param optimizer Pointer to the OptimalTrading instance that is used to
     * find @f[ d* @f] - the optimal entry level.
     * @param model Pointer to the stochastic model to use in the function to
     * evluate @f[ b* @f].
     * @param hitting_time_kernel Pointer to the kernel to use in the function
     * to evaluate @f[ d* @f].
     * @param d_star The upper optimal entry level @f[ d* @f].
     * @param b_star The optimal exit level @f[ b* @f].
     * @param r The discount rate to apply to the optimal mean reversion trading
     * problem.
     * @param c The cost of trading.
     * @return const double The lower bound optimal trading entry level @f[ a*
     * @f].
     */
    const double optimalEntryLower(
        const OptimalTrading* optimizer, const StochasticModel* model,
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& d_star, const double& b_star, const double& r,
        const double& c) const;
    /**
     * @brief Calculates the lower bound optimal entry level @f[ a* @f]
     * for an optimal trading strategy when a stop loss level is provided.
     *
     * @param optimizer Pointer to the OptimalTrading instance that is used to
     * find @f[ d* @f] - the optimal entry level.
     * @param model Pointer to the stochastic model to use in the function to
     * evluate @f[ b* @f].
     * @param hitting_time_kernel Pointer to the kernel to use in the function
     * to evaluate @f[ d* @f].
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
        const OptimalTrading* optimizer, const StochasticModel* model,
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& d_star, const double& b_star, const double& stop_loss,
        const double& r, const double& c) const;
    /**
     * @brief Calculates the optimal entry level @f[ d* @f]
     * for an optimal trading strategy.
     *
     * @param optimizer Pointer to the OptimalTrading instance that is used to
     * find @f[ d* @f] - the optimal entry level.
     * @param model Pointer to the stochastic model to use in the function to
     * evluate @f[ b* @f].
     * @param hitting_time_kernel Pointer to the kernel to use in the function
     * to evaluate @f[ d* @f].
     * @param b_star The optimal exit level @f[ b* @f].
     * @param r The discount rate to apply to the optimal mean reversion trading
     * problem.
     * @param c The cost of trading.
     * @return const double The optimal trading entry level @f[ d* @f].
     */
    const double optimalEntry(
        const OptimalTrading* optimizer, const StochasticModel* model,
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& b_star, const double& r, const double& c) const;
};
#endif  // TRADING_LEVELS_H
