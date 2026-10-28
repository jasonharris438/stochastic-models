#include "stochastic_models/entrypoints/optimal_trading_levels.h"
#include "stochastic_models/exceptions/errors.h"
#include "support/expected_values.h"
#include "support/trading_fixtures.h"

#include <gtest/gtest.h>

namespace {

  using OptimalTradingLevelsTest =
      test_support::TradingFixture<test_support::ornstein_uhlenbeck_set>;
  using OptimalTradingLevelsExponentialTest =
      test_support::TradingFixture<test_support::exponential_set>;
  using OptimalTradingLevelsNoSolutionTest =
      test_support::TradingFixture<test_support::no_solution_set>;

} // namespace

/**
 * @test Tests the output of the optimalEntryLevelLower with a stop loss
 * provided function and asserts that it is near the expected value.
 *
 */
TEST_F(OptimalTradingLevelsTest, optimalEntryLowerStopLossOutputTest) {
  // Declare and initialize test parameters.
  const double stop_loss = 0.04;
  const double d_star = 0.13093;
  const double b_star = 0.455191;
  const double tolerance = 2e-5;

  // Calculate a*.
  const double value = optimalEntryLevelLower(
      d_star, b_star, mu, alpha, sigma, stop_loss, rate, cost
  );

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::optimal_trading_levels_test::
          optimal_entry_lower_stop_loss_output_test,
      tolerance
  ) << "Value produced by optimalEntryLevelLower function "
       "with a stop loss is not equal to the expected value.";
}
/**
 * @test Tests that the optimalEntryLevel function errors due to no solution
 * being found with the provided parameters.
 */
TEST_F(OptimalTradingLevelsNoSolutionTest, optimalEntryLevelNoSolutionTest) {
  // Declare and initialize test parameters.
  const double b_star = 0.750895;

  // Assert that the method is not implemented.
  ASSERT_THROW(
      optimalEntryLevel(b_star, mu, alpha, sigma, rate, cost),
      std::runtime_error
  ) << "optimalEntryLevel did not throw "
       "with invalid parameters.";
}
/**
 * @test Tests the output of the optimalEntryLevel with a stop loss provided
 * function and asserts that it is near the expected value.
 *
 */
TEST_F(OptimalTradingLevelsTest, optimalEntryLevelStopLossOutputTest) {
  // Declare and initialize test parameters.
  const double stop_loss = 0.04;
  const double b_star = 0.455191;
  const double tolerance = 2e-5;

  // Calculate d*.
  const double value =
      optimalEntryLevel(b_star, mu, alpha, sigma, stop_loss, rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::optimal_trading_levels_test::
          optimal_entry_level_stop_loss_output_test,
      tolerance
  ) << "Value produced by optimalEntryLevel function "
       "with a stop loss is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalExitLevel with a stop loss provided
 * function and asserts that it is near the expected value.
 *
 */
TEST_F(OptimalTradingLevelsTest, optimalExitLevelStopLossOutputTest) {
  // Declare and initialize test parameters.
  const double stop_loss = 0.04;
  const double tolerance = 5e-5;

  // Calculate b*.
  const double value =
      optimalExitLevel(mu, alpha, sigma, stop_loss, rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::optimal_trading_levels_test::
          optimal_exit_level_stop_loss_output_test,
      tolerance
  ) << "Value produced by optimalExitLevel function "
       "with a stop loss is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalExitLevel function and asserts
 * that it is near the expected value.
 *
 */
TEST_F(OptimalTradingLevelsTest, optimalExitLevelOutputTest) {
  // Declare and initialize test parameters.
  const double tolerance = 5e-5;

  // Calculate b*.
  const double value = optimalExitLevel(mu, alpha, sigma, rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::optimal_trading_levels_test::optimal_exit_level_output_test,
      tolerance
  ) << "Value produced by optimalExitLevel function "
       "is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalExitLevelExponential function and
 * asserts that it is near the expected value.
 *
 */
TEST_F(
    OptimalTradingLevelsExponentialTest, optimalExitLevelExponentialOutputTest
) {
  // Declare and initialize test parameters.
  const double tolerance = 2e-4;

  // Calculate b*.
  const double value =
      optimalExitLevelExponential(mu, alpha, sigma, rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::optimal_trading_levels_exponential_test::
          optimal_exit_level_exponential_output_test,
      tolerance
  ) << "Value produced by optimalExitLevelExponential function "
       "is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalEntryLevelExponential function and
 * asserts that it is near the expected value.
 *
 */
TEST_F(
    OptimalTradingLevelsExponentialTest, optimalEntryLevelExponentialOutputTest
) {
  // Declare and initialize test parameters.
  const double b_star = 1.4093;
  const double tolerance = 2e-4;

  // Calculate d*.
  const double value =
      optimalEntryLevelExponential(b_star, mu, alpha, sigma, rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::optimal_trading_levels_exponential_test::
          optimal_entry_level_exponential_output_test,
      tolerance
  ) << "Value produced by optimalEntryLevelExponential function "
       "is not equal to the expected value.";
}
/**
 * @test Tests that the optimalExitLevelExponential function errors due to
 * no solution being found with the provided parameters.
 *
 */
TEST_F(
    OptimalTradingLevelsNoSolutionTest,
    optimalExitLevelExponentialNoSolutionTest
) {
  // Assert that the method is not implemented.
  ASSERT_THROW(
      optimalExitLevelExponential(mu, alpha, sigma, rate, cost),
      std::runtime_error
  ) << "optimalExitLevelExponential did not throw a runtime error when "
       "provided "
       "with invalid parameters.";
}
/**
 * @test Tests that the optimalEntryLevelLowerExponential function throws
 * RootNotBracketedError when the d*-equation has no root in the bracket.
 *
 */
TEST_F(
    OptimalTradingLevelsExponentialTest,
    optimalEntryLevelLowerExponentialNoRootTest
) {
  // Declare and initialize test parameters.
  const double d_star = 1.236887;
  const double b_star = 1.4093;

  // Assert that the no-root bracket is reported.
  EXPECT_THROW(
      optimalEntryLevelLowerExponential(
          d_star, b_star, mu, alpha, sigma, rate, cost
      ),
      RootNotBracketedError
  ) << "optimalEntryLevelLowerExponential did not throw "
       "RootNotBracketedError for a bracket with no root.";
}
/**
 * @test Tests the output of the optimalEntryLevel function and asserts
 * that it is near the expected value.
 *
 */
TEST_F(OptimalTradingLevelsTest, optimalEntryLevelOutputTest) {
  // Declare and initialize test parameters.
  const double b_star = 0.466836;
  const double tolerance = 2e-5;

  // Calculate d*.
  const double value = optimalEntryLevel(b_star, mu, alpha, sigma, rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::optimal_trading_levels_test::optimal_entry_level_output_test,
      tolerance
  ) << "Value produced by optimalEntryLevel function "
       "is not equal to the expected value.";
}
