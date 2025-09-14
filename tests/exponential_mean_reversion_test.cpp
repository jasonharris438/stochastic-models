#include "stochastic_models/trading/exponential_mean_reversion.h"

#include <gtest/gtest.h>

#include <cstdlib>
#include <stdexcept>

#include "stochastic_models/sde/ornstein_uhlenbeck.h"
#include "stochastic_models/trading/optimal_mean_reversion.h"
/**
 * @test Tests the method ExponentialMeanReversion::b with a stop_loss argument
 * is not implemented
 *
 */
TEST(ExponentialMeanReversionTest, stopLossBNotImplementedTest) {
    // Declare and initialize model and test parameters.
    const double x = 0.2;
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double stop_loss = 0.05;
    const double c = 0.02;
    const double r = 0.05;

    // Create core model and optimal mean reversion instances.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
    ExponentialMeanReversion mean_reversion;

    // Assert that the method is not implemented.
    ASSERT_THROW(mean_reversion.b(x, hitting_time_kernel, stop_loss, r, c),
                 std::logic_error)
        << "OptimalMeanReversion::b with stop loss is not meant to be "
           "implemented.";
    delete hitting_time_kernel;
}
/**
 * @test Tests the method ExponentialMeanReversion::d with a stop_loss argument
 * is not implemented
 *
 */
TEST(ExponentialMeanReversionTest, stopLossDNotImplementedTest) {
    // Declare and initialize model and test parameters.
    const double x = 0.2;
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double b_star = 0.466;
    const double stop_loss = 0.05;
    const double c = 0.02;
    const double r = 0.05;

    // Create core model and optimal mean reversion instances.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
    ExponentialMeanReversion mean_reversion;

    // Assert that the method is not implemented.
    ASSERT_THROW(
        mean_reversion.d(x, hitting_time_kernel, b_star, stop_loss, r, c),
        std::logic_error)
        << "OptimalMeanReversion::d with stop loss is not meant to be "
           "implemented.";
    delete hitting_time_kernel;
}
/**
 * @test Tests the method ExponentialMeanReversion::a with a stop_loss argument
 * is not implemented
 *
 */
TEST(ExponentialMeanReversionTest, stopLossANotImplementedTest) {
    // Declare and initialize model and test parameters.
    const double x = 0.2;
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double b_star = 0.466;
    const double stop_loss = 0.05;
    const double c = 0.02;
    const double r = 0.05;

    // Create core model and optimal mean reversion instances.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
    ExponentialMeanReversion mean_reversion;

    // Assert that the method is not implemented.
    ASSERT_THROW(
        mean_reversion.a(x, hitting_time_kernel, b_star, stop_loss, r, c),
        std::logic_error)
        << "OptimalMeanReversion::a with stop loss is not meant to be "
           "implemented.";
    delete hitting_time_kernel;
}
/**
 * @test Tests the method ExponentialMeanReversion::V with a stop_loss argument
 * is not implemented
 *
 */
TEST(ExponentialMeanReversionTest, stopLossVNotImplementedTest) {
    // Declare and initialize model and test parameters.
    const double x = 0.2;
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double b_star = 0.466;
    const double stop_loss = 0.05;
    const double c = 0.02;
    const double r = 0.05;
    // Create core model and optimal mean reversion instances.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
    ExponentialMeanReversion mean_reversion;

    // Assert that the method is not implemented.
    ASSERT_THROW(
        mean_reversion.V(hitting_time_kernel, x, b_star, stop_loss, r, c),
        std::logic_error)
        << "OptimalMeanReversion::V with stop loss is not meant to be "
           "implemented.";
    delete hitting_time_kernel;
}
/**
 * @test Tests the method ExponentialMeanReversion::instantaneousDifferential
 * with a stop_loss argument is not implemented
 *
 */
TEST(ExponentialMeanReversionTest,
     stopLossinstantaneousDifferentialNotImplementedTest) {
    // Declare and initialize model and test parameters.
    const double x = 0.2;
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double b_star = 0.466;
    const double stop_loss = 0.05;
    const double c = 0.02;
    const double r = 0.05;
    ModelFunc fn = funcIntegrateF;

    // Create core model and optimal mean reversion instances.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
    ExponentialMeanReversion mean_reversion;

    // Assert that the method is not implemented.
    ASSERT_THROW(mean_reversion.instantaneousDifferential(
                     fn, hitting_time_kernel, x, b_star, stop_loss, r, c),
                 std::logic_error)
        << "OptimalMeanReversion::instantaneousDifferential with stop loss is "
           "not meant to be implemented.";
    delete hitting_time_kernel;
}
