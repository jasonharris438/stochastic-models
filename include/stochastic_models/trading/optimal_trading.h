#ifndef OPTIMAL_TRADING_H
#define OPTIMAL_TRADING_H
#include "stochastic_models/hitting_times/hitting_time_ornstein_uhlenbeck.h"
#include "stochastic_models/numeric_utils/types.h"
#include "stochastic_models/sde/stochastic_model.h"

/**
 * @brief The OptimalTrading class is an abstract class that defines methods
 * that are used in calculating the optimal trading strategy.
 *
 */
class OptimalTrading {
   public:
    /**
     * @brief Construct a new OptimalTrading object and return on heap memory
     * using the class' copy constructor in the caller instance
     *
     * @return const OptimalTrading* Pointer to the new instance.
     */
    virtual const OptimalTrading* clone() const = 0;
    /**
     * @brief Abstract method that calculates @f[ F(x;r) @f] for the optimal
     * trading strategy.
     *
     * @param model Pointer to the stochastic model to use in the function @f[
     * F(x;r) @f].
     * @param x The current value x to evaluate @f[ F(x;r) @f].
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value of @f[ F(x;r) @f].
     */
    virtual const double F(
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& x, const double& r, const double& c) const = 0;
    /**
     * @brief Abstract method that calculates @f[ G(x;r) @f] for the optimal
     * trading strategy.
     *
     * @param model Pointer to the stochastic model to use in the function @f[
     * G(x;r) @f].
     * @param x The current value x to evaluate @f[ G(x;r) @f].
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value of @f[ G(x;r) @f].
     */
    virtual const double G(
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& x, const double& r, const double& c) const = 0;
    /**
     * @brief Abstract method that calculates @f[ L^* @f] for the optimal
     * trading strategy.
     *
     * @param model Pointer to the stochastic model to use in the function @f[
     * L^* @f].
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value of @f[ L^* @f].
     */
    virtual const double L_star(
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& r, const double& c) const = 0;
    /**
     * @brief Abstract method that calculates the value @f[ b @f]
     * for the optimal trading strategy.
     *
     * @param value The value at which to evaluate the optimal exit level.
     * @param model Pointer to the stochastic model to use in the function @f[ b
     * @f].
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value @f[ b @f].
     */
    virtual const double b(
        const double& value,
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& r, const double& c) const = 0;
    /**
     * @brief Abstract method that calculates the value @f[ b @f]
     * for the optimal trading strategy with a stop loss level.
     *
     * @param value The value at which to evaluate the optimal exit level.
     * @param model Pointer to the stochastic model to use in the function @f[ b
     * @f].
     * @param stop_loss The stop loss level to use in the optimal trading.
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value @f[ b @f].
     */
    virtual const double b(
        const double& value,
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& stop_loss, const double& r, const double& c) const = 0;
    /**
     * @brief Abstract method that calculates the value @f[ a @f], the lower
     * bound for the optimal entry region for the optimal trading strategy.
     *
     * @param value The value at which to evaluate the optimal entry level.
     * @param model Pointer to the stochastic model to use in the function @f[ d
     * @f].
     * @param b_star The optimal exit level calculated by the OptimalTrading
     * instance.
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value @f[ a @f].
     */
    virtual const double a(
        const double& value,
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& b_star, const double& r, const double& c) const = 0;
    /**
     * @brief Abstract method that calculates the value @f[ a @f], the lower
     * bound for the optimal entry region for the optimal trading strategy when
     * a stop loss level is provided.
     *
     * @param value The value at which to evaluate the optimal entry level.
     * @param model Pointer to the stochastic model to use in the function @f[ d
     * @f].
     * @param b_star The optimal exit level calculated by the OptimalTrading
     * instance.
     * @param stop_loss The stop loss level to use in the optimal trading.
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value @f[ a @f].
     */
    virtual const double a(
        const double& value,
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& b_star, const double& stop_loss, const double& r,
        const double& c) const = 0;
    /**
     * @brief Abstract method that calculates the value @f[ d @f] for the
     * optimal trading strategy when a stop loss level is provided.
     *
     * @param value The value at which to evaluate the optimal entry level.
     * @param model Pointer to the stochastic model to use in the function @f[ d
     * @f].
     * @param b_star The optimal exit level calculated by the OptimalTrading
     * instance.
     * @param stop_loss The stop loss level to use in the optimal trading.
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value @f[ d @f].
     */
    virtual const double d(
        const double& value,
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& b_star, const double& stop_loss, const double& r,
        const double& c) const = 0;
    /**
     * @brief Abstract method that calculates the value @f[ d @f] for the
     * optimal trading strategy.
     *
     * @param value The value at which to evaluate the optimal entry level.
     * @param model Pointer to the stochastic model to use in the function @f[ d
     * @f].
     * @param b_star The optimal exit level calculated by the OptimalTrading
     * instance.
     * @param r The discount rate to apply to the optimal trading problem.
     * @param c The cost of trading.
     * @return const double The value @f[ d @f].
     */
    virtual const double d(
        const double& value,
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& b_star, const double& r, const double& c) const = 0;
    /**
     * @brief Abstract method that calculates the instantaneous rate of change
     * of the function fn at point x using the model pointed to by model. This
     * overloaded version requires a stop loss level to be provided.
     *
     * @param fn The function with which to calculate the instantaneous rate of
     * change.
     * @param model The model to use in calculating the instantaneous rate of
     * change.
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
    virtual const double instantaneousDifferential(
        ModelFunc fn, const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& x, const double& b_star, const double& stop_loss,
        const double& r, const double& c) const = 0;
    /**
     * @brief Abstract method that calculates the instantaneous rate of change
     * of the function fn at point x using the model pointed to by model.
     *
     * @param fn The function with which to calculate the instantaneous rate of
     * change.
     * @param model The model to use in calculating the instantaneous rate of
     * change.
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
    virtual const double instantaneousDifferential(
        ModelFunc fn, const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& x, const double& b_star, const double& r,
        const double& c) const = 0;
    /**
     * @brief Abstract method that calculates the instantaneous rate of change
     * of the function fn at point x using the model pointed to by model.
     *
     * @param fn The function with which to calculate the instantaneous rate of
     * change.
     * @param model The model to use in calculating the instantaneous rate of
     * change.
     * @param x The point at which to calculate the instantaneous rate of
     * change.
     * @param r The discount rate to apply to provide to the function when
     * differentiating.
     * @param c The cost of trading to provide to the function when
     * differentiating.
     * @return const double The instantaneous rate of change of the function at
     * the point x.
     */
    virtual const double instantaneousDifferential(
        ModelFunc fn, const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& x, const double& r, const double& c) const = 0;
    /**
     * @brief Abstract method that determines the value function for the optimal
     * trading model.
     *
     * @param model The model to use in calculating the value function.
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
    virtual const double V(
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& x, const double& b_star, const double& r,
        const double& c) const = 0;
    /**
     * @brief Abstract method that determines the value function for the optimal
     * trading model when a stop loss is provided.
     *
     * @param model The model to use in calculating the value function.
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
    virtual const double V(
        const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
        const double& x, const double& b_star, const double& stop_loss,
        const double& r, const double& c) const = 0;
};
#endif  // OPTIMAL_TRADING_H
