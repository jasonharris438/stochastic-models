#ifndef OPTIMAL_MEAN_REVERSION_H
#define OPTIMAL_MEAN_REVERSION_H
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
 * value function @f[ V(x) @f].
 *
 * @param x The point at which to evaluate the function @f[ V(x) @f].
 * @param params The parameters provided to use to implement the function @f[
 * V(x) @f] with respect to the model type.
 * @return double The value of the integral @f[ V(x) @f] at the point x.
 */
double valueFunc(double x, void* params);
/**
 * @brief Used when evaluating the optimal mean reversion trading
 * value function @f[ V(x) @f] when a stop loss is provided.
 *
 * @param x The point at which to evaluate the function @f[ V(x) @f].
 * @param params The parameters provided to use to implement the function @f[
 * V(x) @f] with respect to the model type.
 * @return double The value of the integral @f[ V(x) @f] at the point x.
 */
double valueFuncStopLoss(double x, void* params);
/**
 * @brief Used when integrating to evaluate the optimal mean reversion trading
 * model function @f[ G(x;r) @f].
 *
 * @param x The point at which to evaluate the function @f[ G(x;r) @f].
 * @param params The parameters provided to use to implement the function @f[
 * G(x;r) @f] with respect to the model type.
 * @return double The value of the integral @f[ G(x;r) @f] at the point x.
 */
double funcIntegrateG(double x, void* params);
/**
 * @brief Used when integrating to evaluate the optimal mean reversion trading
 * model function @f[ F(x;r) @f].
 *
 * @param x The point at which to evaluate the function @f[ F(x;r) @f].
 * @param params The parameters provided to use to implement the function @f[
 * F(x;r) @f] with respect to the model type.
 * @return double The value of the integral @f[ F(x;r) @f] at the point x.
 */
double funcIntegrateF(double x, void* params);
/**
 * @brief Used when calculating the optimal mean reversion trading model
 * function @f[ F(x;r) @f].
 *
 * @param x The point at which to evaluate the antiderivative.
 * @param params The parameters provided via the numerical integration function
 * that contains the model and discount rate to use.
 * @return double The value of the integral @f[ F(x;r) @f] at the point x.
 */
double funcOptimalMeanReversionF(double x, void* params);
/**
 * @brief Used when calculating the optimal mean reversion trading model
 * function @f[ G(x;r) @f].
 *
 * @param x The point at which to evaluate the antiderivative.
 * @param params The parameters provided via the numerical integration function
 * that contains the model and discount rate to use.
 * @return double The value of the integral @f[ G(x;r) @f] at the point x.
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
     * @brief Calculates @f[ F(x;r) @f] in the optimal
     * trading strategy for a mean reverting model.
     *
     * @param hitting_time_kernel Pointer to the hitting time kernel instance to
     * use in the function @f[ F(x;r) @f].
     * @param x The current value x to evaluate @f[ F(x;r) @f].
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value of @f[ F(x;r) @f].
     */
    const double F(const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
                   const double& x, const double& r,
                   const double& c) const override;
    /**
     * @brief Calculates @f[ G(x;r) @f] in the optimal
     * trading strategy for a mean reverting model.
     *
     * @param hitting_time_kernel Pointer to the hitting time kernel instance to
     * use in the function @f[ G(x;r) @f].
     * @param x The current value x to evaluate @f[ G(x;r) @f].
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value of @f[ G(x;r) @f].
     */
    const double G(const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
                   const double& x, const double& r,
                   const double& c) const override;
    /**
     * @brief Calculates the value @f[ L^* @f] for the optimal trading strategy.
     *
     * @param hitting_time_kernel Pointer to the hitting time kernel to use in
     * the function
     * @f[ L^* @f].
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value of @f[ L^* @f].
     */
    const double L_star(const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
                        const double& r, const double& c) const override;
    /**
     * @brief Calculates the value @f[ b @f] for the optimal trading strategy
     * for a mean reverting model with a stop loss level.
     *
     * @param value The value at which to evaluate the optimal exit level.
     * @param hitting_time_kernel Pointer to the hitting time kernel to use in
     * the function @f[ b
     * @f].
     * @param stop_loss The stop loss level to use in the optimal trading.
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value @f[ b @f].
     */
    const double b(const double& value,
                   const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
                   const double& stop_loss, const double& r,
                   const double& c) const override;
    /**
     * @brief Calculates the value @f[ b @f] in the optimal trading strategy for
     * a mean reverting model.
     *
     * @param value The value at which to evaluate the optimal exit level.
     * @param hitting_time_kernel Pointer to the hitting time kernel to use in
     * the function @f[ b
     * @f].
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value @f[ b @f].
     */
    const double b(const double& value,
                   const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
                   const double& r, const double& c) const override;
    /**
     * @brief Calculates the value @f[ a @f], the lower bound for the optimal
     * entry region for the optimal trading strategy for a mean reverting model.
     *
     * @param value The value at which to evaluate the optimal entry level.
     * @param hitting_time_kernel Pointer to the hitting time kernel to use in
     * the function @f[ d
     * @f].
     * @param b_star The optimal exit level calculated by the OptimalTrading
     * instance.
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value @f[ a @f].
     */
    const double a(const double& value,
                   const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
                   const double& b_star, const double& r,
                   const double& c) const override;
    /**
     * @brief Calculates the value @f[ a @f], the lower bound for the optimal
     * entry region for the optimal trading strategy for a mean reverting model
     * when a stop loss level is provided.
     *
     * @param value The value at which to evaluate the optimal entry level.
     * @param hitting_time_kernel Pointer to the hitting time kernel to use in
     * the function @f[ d
     * @f].
     * @param b_star The optimal exit level calculated by the OptimalTrading
     * instance.
     * @param stop_loss The stop loss level to use in the optimal trading.
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value @f[ a @f].
     */
    const double a(const double& value,
                   const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
                   const double& b_star, const double& stop_loss,
                   const double& r, const double& c) const override;
    /**
     * @brief Calculates the value @f[ d @f] for the optimal trading strategy
     * for a mean reverting model when a stop loss level is provided.
     *
     * @param value The value at which to evaluate the optimal entry level.
     * @param hitting_time_kernel Pointer to the hitting time kernel to use in
     * the function @f[ d
     * @f].
     * @param b_star The optimal exit level calculated by the OptimalTrading
     * instance.
     * @param stop_loss The stop loss level to use in the optimal trading.
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value @f[ d @f].
     */
    const double d(const double& value,
                   const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
                   const double& b_star, const double& stop_loss,
                   const double& r, const double& c) const override;
    /**
     * @brief Calculates the value @f[ d @f] for the optimal trading strategy
     * for a mean reverting model.
     *
     * @param value The value at which to evaluate the optimal entry level.
     * @param hitting_time_kernel Pointer to the hitting time kernel to use in
     * the function @f[ d
     * @f].
     * @param b_star The optimal exit level calculated by the OptimalTrading
     * instance.
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value @f[ d @f].
     */
    const double d(const double& value,
                   const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
                   const double& b_star, const double& r,
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
        ModelFunc fn, const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& x, const double& r, const double& c) const override;
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
        ModelFunc fn, const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& x, const double& b_star, const double& r,
        const double& c) const override;
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
        ModelFunc fn, const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& x, const double& b_star, const double& stop_loss,
        const double& r, const double& c) const override;
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
    const double V(const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
                   const double& x, const double& b_star, const double& r,
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
    const double V(const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
                   const double& x, const double& b_star,
                   const double& stop_loss, const double& r,
                   const double& c) const override;
};
#endif  // OPTIMAL_MEAN_REVERSION_H
