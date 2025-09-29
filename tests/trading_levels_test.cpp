#include "stochastic_models/trading/trading_levels.h"

#include <gtest/gtest.h>

#include <cstdlib>

#include "stochastic_models/numeric_utils/helpers.h"
#include "stochastic_models/numeric_utils/integration.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"
#include "stochastic_models/trading/exponential_mean_reversion.h"
#include "stochastic_models/trading/optimal_mean_reversion.h"
#include "stochastic_models/trading/trading_levels_exponential.h"
#include "stochastic_models/trading/trading_levels_params.h"
/**
 * @test Tests the output of the TradingLevels::optimalExit method with a stop
 * loss provided and asserts that it is near the expected value.
 *
 */
TEST(TradingLevelsTest, exitLevelStopLossOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double stop_loss = 0.05;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-5;

    // Create trading levels instance to manage allocations.
    OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

    // Calculate b*.
    const double value = tradingLevels.optimalExit(stop_loss, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.46912012), tolerance)
        << "Value produced by OrnsteinUhlenbeckTradingLevels::optimalExit "
           "with a stop loss provided is not equal to the expected value.";
}
/**
 * @test Tests the output of the TradingLevels::optimalExit method and asserts
 * that it is near the expected value.
 *
 */
TEST(TradingLevelsTest, exitLevelOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-5;

    // Create trading levels instance to manage allocations.
    OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

    // Calculate b*.
    const double value = tradingLevels.optimalExit(r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.466836), tolerance)
        << "Value produced by OrnsteinUhlenbeckTradingLevels::optimalExit "
           "is not equal to the expected value.";
}
/**
 * @test Tests the output of the TradingLevels::optimalExit with
 * ExponentialMeanReversion optimizer instance method and asserts that it is
 * near the expected value.
 *
 */
TEST(TradingLevelsTest, exitLevelExponentialOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 5;
    const double mu = 1.3499;
    const double sigma = 0.15;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Create trading levels instance to manage allocations.
    OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

    // Calculate b*.
    const double value = tradingLevels.optimalExit(r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 1.4093), tolerance)
        << "Value produced by "
           "OrnsteinUhlenbeckTradingLevelsExponential::optimalExit "
           "with ExponentialMeanReversion optimizer is not equal to the "
           "expected value.";
}
/**
 * @test Tests the output of the TradingLevels::optimalEntryLower method and
 * asserts that it is near the expected value when a stop loss level is
 * provided.
 *
 */
TEST(TradingLevelsTest, entryLevelLowerBoundStopLossOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double stop_loss = 0.05;
    const double d_star = 0.136755;
    const double b_star = 0.450895;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Create trading levels instance to manage allocations.
    OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

    // Calculate d*.
    const double value =
        tradingLevels.optimalEntryLower(d_star, b_star, stop_loss, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.118451), tolerance)
        << "Value produced by "
           "OrnsteinUhlenbeckTradingLevels::optimalEntryLower "
           "is not equal to the expected value when a stop loss is provided.";
}
/**
 * @test Tests the output of the TradingLevels::optimalEntry method and asserts
 * that it is near the expected value when a stop loss level is provided.
 *
 */
TEST(TradingLevelsTest, entryLevelStopLossOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double stop_loss = 0.05;
    const double b_star = 0.450895;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Create trading levels instance to manage allocations.
    OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

    // Calculate d*.
    const double value = tradingLevels.optimalEntry(b_star, stop_loss, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.136755), tolerance)
        << "Value produced by OrnsteinUhlenbeckTradingLevels::optimalEntry "
           "is not equal to the expected value when a stop loss is provided.";
}
/**
 * @test Tests the output of the TradingLevels::optimalEntry method and asserts
 * that it is near the expected value.
 *
 */
TEST(TradingLevelsTest, entryLevelOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double b_star = 0.466836;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Create trading levels instance to manage allocations.
    OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

    // Calculate d*.
    const double value = tradingLevels.optimalEntry(b_star, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.116948), tolerance)
        << "Value produced by OrnsteinUhlenbeckTradingLevels::optimalEntry "
           "is not equal to the expected value.";
}
/**
 * @test Tests the output of the TradingLevels::optimalEntry method and asserts
 * that it throws safely when the upper and lower bound values are invalid.
 *
 */
TEST(TradingLevelsTest, entryLevelBoundsErrorTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double b_star = -4.466836;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Create trading levels instance to manage allocations.
    OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

    // Assert that the method is not implemented.
    ASSERT_THROW(tradingLevels.optimalEntry(b_star, r, c),
                 std::invalid_argument)
        << "OrnsteinUhlenbeckTradingLevels::optimalEntry did not throw "
           "with invalid upper and lower bound values.";
}
/**
 * @test Tests the output of the TradingLevels::optimalEntry method with an
 * ExponentialMeanReversion optimizer instance and asserts that it is near the
 * expected value.
 *
 */
TEST(TradingLevelsTest, entryLevelExponentialOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 5;
    const double mu = 1.3499;
    const double sigma = 0.15;
    const double b_star = 1.4093;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Create trading levels instance to manage allocations.
    OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

    // Calculate d*.
    const double value = tradingLevels.optimalEntry(b_star, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 1.24096), tolerance)
        << "Value produced by "
           "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntry "
           "with ExponentialMeanReversion optimizer is not equal to the "
           "expected value.";
}
/**
 * @test Tests the output of the TradingLevels::optimalEntryLower method with an
 * ExponentialMeanReversion optimizer instance and asserts that it is near the
 * expected value.
 *
 */
TEST(TradingLevelsTest, entryLevelLowerExponentialOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 5;
    const double mu = 1.3499;
    const double sigma = 0.15;
    const double d_star = 1.24096;
    const double b_star = 1.4093;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Create trading levels instance to manage allocations.
    OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

    // Calculate d*.
    const double value = tradingLevels.optimalEntryLower(d_star, b_star, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 1.16016), tolerance)
        << "Value produced by "
           "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntryLower "
           "with ExponentialMeanReversion optimizer is not equal to the "
           "expected value.";
}
