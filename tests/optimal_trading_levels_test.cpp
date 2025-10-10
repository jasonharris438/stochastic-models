#include "stochastic_models/entrypoints/optimal_trading_levels.h"

#include <gtest/gtest.h>

#include <cstdlib>

#include "stochastic_models/numeric_utils/helpers.h"

/**
 * @test Tests the output of the optimalEntryLevelLower with a stop loss
 * provided function and asserts that it is near the expected value.
 *
 */
TEST(OptimalTradingLevelsTest, optimalEntryLowerStopLossOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double stop_loss = 0.05;
    const double d_star = 0.136755;
    const double b_star = 0.450895;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-5;

    // Calculate a*.
    const double value = optimalEntryLevelLower(d_star, b_star, mu, alpha,
                                                sigma, stop_loss, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.118451), tolerance)
        << "Value produced by optimalEntryLevelLower function "
           "with a stop loss is not equal to the expected value.";
}
/**
 * @test Tests that the optimalEntryLevel function errors due to no solution
 * being found with the provided parameters.
 */
TEST(OptimalTradingLevelsTest, optimalEntryLevelNoSolutionTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 0.000116;
    const double mu = 1.818978;
    const double sigma = 0.006623;
    const double b_star = 0.750895;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-5;

    // Assert that the method is not implemented.
    ASSERT_THROW(optimalEntryLevel(b_star, mu, alpha, sigma, r, c),
                 std::runtime_error)
        << "optimalEntryLevel did not throw "
           "with invalid parameters.";
}
/**
 * @test Tests the output of the optimalEntryLevel with a stop loss provided
 * function and asserts that it is near the expected value.
 *
 */
TEST(OptimalTradingLevelsTest, optimalEntryLevelStopLossOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double stop_loss = 0.05;
    const double b_star = 0.450895;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-5;

    // Calculate d*.
    const double value =
        optimalEntryLevel(b_star, mu, alpha, sigma, stop_loss, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.136755), tolerance)
        << "Value produced by optimalEntryLevel function "
           "with a stop loss is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalExitLevel with a stop loss provided
 * function and asserts that it is near the expected value.
 *
 */
TEST(OptimalTradingLevelsTest, optimalExitLevelStopLossOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double stop_loss = 0.05;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-5;

    // Calculate b*.
    const double value = optimalExitLevel(mu, alpha, sigma, stop_loss, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.46912012), tolerance)
        << "Value produced by optimalExitLevel function "
           "with a stop loss is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalExitLevel function and asserts
 * that it is near the expected value.
 *
 */
TEST(OptimalTradingLevelsTest, optimalExitLevelOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-5;

    // Calculate b*.
    const double value = optimalExitLevel(mu, alpha, sigma, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.466836), tolerance)
        << "Value produced by optimalExitLevel function "
           "is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalExitLevelExponential function and
 * asserts that it is near the expected value.
 *
 */
TEST(OptimalTradingLevelsTest, optimalExitLevelExponentialOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 5;
    const double mu = 1.3499;
    const double sigma = 0.15;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Calculate b*.
    const double value = optimalExitLevelExponential(mu, alpha, sigma, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 1.4093), tolerance)
        << "Value produced by optimalExitLevelExponential function "
           "is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalEntryLevelExponential function and
 * asserts that it is near the expected value.
 *
 */
TEST(OptimalTradingLevelsTest, optimalEntryLevelExponentialOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 5;
    const double mu = 1.3499;
    const double sigma = 0.15;
    const double b_star = 1.4093;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Calculate d*.
    const double value =
        optimalEntryLevelExponential(b_star, mu, alpha, sigma, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 1.24096), tolerance)
        << "Value produced by optimalEntryLevelExponential function "
           "is not equal to the expected value.";
}
/**
 * @test Tests that the optimalExitLevelExponential function errors due to
 * no solution being found with the provided parameters.
 *
 */
TEST(OptimalTradingLevelsTest, optimalExitLevelExponentialNoSolutionTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 0.000116;
    const double mu = 1.818978;
    const double sigma = 0.006623;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Assert that the method is not implemented.
    ASSERT_THROW(optimalExitLevelExponential(mu, alpha, sigma, r, c),
                 std::runtime_error)
        << "optimalExitLevelExponential did not throw a runtime error when "
           "provided "
           "with invalid parameters.";
}
/**
 * @test Tests the output of the optimalEntryLevelLowerExponential function and
 * asserts that it is near the expected value.
 *
 */
TEST(OptimalTradingLevelsTest, optimalEntryLevelLowerExponentialOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 5;
    const double mu = 1.3499;
    const double sigma = 0.15;
    const double d_star = 1.24096;
    const double b_star = 1.4093;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Calculate d*.
    const double value = optimalEntryLevelLowerExponential(d_star, b_star, mu,
                                                           alpha, sigma, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 1.16016), tolerance)
        << "Value produced by optimalEntryLevelLowerExponential function "
           "is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalEntryLevel function and asserts
 * that it is near the expected value.
 *
 */
TEST(OptimalTradingLevelsTest, optimalEntryLevelOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double b_star = 0.466836;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Calculate d*.
    const double value = optimalEntryLevel(b_star, mu, alpha, sigma, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.116948), tolerance)
        << "Value produced by optimalEntryLevel function "
           "is not equal to the expected value.";
}
