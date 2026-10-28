#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/numeric_utils/integration.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"
#include "stochastic_models/trading/exponential_mean_reversion.h"
#include "stochastic_models/trading/optimal_mean_reversion.h"
#include "stochastic_models/trading/trading_levels.h"
#include "stochastic_models/trading/trading_levels_exponential.h"
#include "stochastic_models/trading/trading_levels_params.h"
#include "support/expected_values.h"
#include "support/trading_fixtures.h"

#include <gtest/gtest.h>
#include <stdexcept>

namespace {

  class TradingLevelsTest : public test_support::TradingFixture<
                                test_support::ornstein_uhlenbeck_set> {
  protected:
    OrnsteinUhlenbeckTradingLevels trading_levels{mu, alpha, sigma};
  };

  class TradingLevelsExponentialTest
      : public test_support::TradingFixture<test_support::exponential_set> {
  protected:
    OrnsteinUhlenbeckTradingLevelsExponential trading_levels{mu, alpha, sigma};
  };

} // namespace

/**
 * @test Tests the output of the TradingLevels::optimalExit method with a stop
 * loss provided and asserts that it is near the expected value.
 *
 */
TEST_F(TradingLevelsTest, exitLevelStopLossOutputTest) {
  // Declare and initialize test parameters.
  const double stop_loss = 0.04;
  const double tolerance = 5e-5;

  // Calculate b*.
  const double value = trading_levels.optimalExit(stop_loss, rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value, expected::trading_levels_test::exit_level_stop_loss_output_test,
      tolerance
  ) << "Value produced by OrnsteinUhlenbeckTradingLevels::optimalExit "
       "with a stop loss provided is not equal to the expected value.";
}
/**
 * @test Tests the output of the TradingLevels::optimalExit method and asserts
 * that it is near the expected value.
 *
 */
TEST_F(TradingLevelsTest, exitLevelOutputTest) {
  // Declare and initialize test parameters.
  const double tolerance = 5e-5;

  // Calculate b*.
  const double value = trading_levels.optimalExit(rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value, expected::trading_levels_test::exit_level_output_test, tolerance
  ) << "Value produced by OrnsteinUhlenbeckTradingLevels::optimalExit "
       "is not equal to the expected value.";
}
/**
 * @test Tests the output of the TradingLevels::optimalExit with
 * ExponentialMeanReversion optimizer instance method and asserts that it is
 * near the expected value.
 *
 */
TEST_F(TradingLevelsExponentialTest, exitLevelExponentialOutputTest) {
  // Declare and initialize test parameters.
  const double tolerance = 2e-4;

  // Calculate b*.
  const double value = trading_levels.optimalExit(rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::trading_levels_exponential_test::
          exit_level_exponential_output_test,
      tolerance
  ) << "Value produced by "
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
TEST_F(TradingLevelsTest, entryLevelLowerBoundStopLossOutputTest) {
  // Declare and initialize test parameters.
  const double stop_loss = 0.04;
  const double d_star = 0.13093;
  const double b_star = 0.455191;
  const double tolerance = 2e-5;

  // Calculate d*.
  const double value =
      trading_levels.optimalEntryLower(d_star, b_star, stop_loss, rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::trading_levels_test::
          entry_level_lower_bound_stop_loss_output_test,
      tolerance
  ) << "Value produced by "
       "OrnsteinUhlenbeckTradingLevels::optimalEntryLower "
       "is not equal to the expected value when a stop loss is provided.";
}
/**
 * @test Tests the output of the TradingLevels::optimalEntry method and asserts
 * that it is near the expected value when a stop loss level is provided.
 *
 */
TEST_F(TradingLevelsTest, entryLevelStopLossOutputTest) {
  // Declare and initialize test parameters.
  const double stop_loss = 0.04;
  const double b_star = 0.455191;
  const double tolerance = 2e-5;

  // Calculate d*.
  const double value =
      trading_levels.optimalEntry(b_star, stop_loss, rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value, expected::trading_levels_test::entry_level_stop_loss_output_test,
      tolerance
  ) << "Value produced by OrnsteinUhlenbeckTradingLevels::optimalEntry "
       "is not equal to the expected value when a stop loss is provided.";
}
/**
 * @test Tests the output of the TradingLevels::optimalEntry method and asserts
 * that it is near the expected value.
 *
 */
TEST_F(TradingLevelsTest, entryLevelOutputTest) {
  // Declare and initialize test parameters.
  const double b_star = 0.466836;
  const double tolerance = 2e-5;

  // Calculate d*.
  const double value = trading_levels.optimalEntry(b_star, rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value, expected::trading_levels_test::entry_level_output_test, tolerance
  ) << "Value produced by OrnsteinUhlenbeckTradingLevels::optimalEntry "
       "is not equal to the expected value.";
}
/**
 * @test Tests the output of the TradingLevels::optimalEntry method and asserts
 * that it throws safely when the upper and lower bound values are invalid.
 *
 */
TEST_F(TradingLevelsTest, entryLevelBoundsErrorTest) {
  // Declare and initialize test parameters.
  const double b_star = -4.466836;

  // Assert that the method is not implemented.
  ASSERT_THROW(
      trading_levels.optimalEntry(b_star, rate, cost), std::invalid_argument
  ) << "OrnsteinUhlenbeckTradingLevels::optimalEntry did not throw "
       "with invalid upper and lower bound values.";
}
/**
 * @test Asserts that the OU trading levels reject a no-stop-loss lower entry
 * level a*, which has no mathematical definition for the Ornstein-Uhlenbeck
 * model.
 *
 */
TEST_F(TradingLevelsTest, entryLevelLowerNotImplementedTest) {
  // Declare and initialize test parameters.
  const double b_star = 0.4551908;
  const double d_star = 0.13093;

  // Assert that the method is not implemented.
  EXPECT_THROW(
      trading_levels.optimalEntryLower(d_star, b_star, rate, cost),
      std::logic_error
  ) << "OrnsteinUhlenbeckTradingLevels::optimalEntryLower without a stop "
       "loss did not throw.";
}
/**
 * @test Tests the output of the TradingLevels::optimalEntry method with an
 * ExponentialMeanReversion optimizer instance and asserts that it is near the
 * expected value.
 *
 */
TEST_F(TradingLevelsExponentialTest, entryLevelExponentialOutputTest) {
  // Declare and initialize test parameters.
  const double b_star = 1.4093;
  const double tolerance = 2e-4;

  // Calculate d*.
  const double value = trading_levels.optimalEntry(b_star, rate, cost);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::trading_levels_exponential_test::
          entry_level_exponential_output_test,
      tolerance
  ) << "Value produced by "
       "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntry "
       "with ExponentialMeanReversion optimizer is not equal to the "
       "expected value.";
}
/**
 * @test Tests that the TradingLevels::optimalEntryLower method with an
 * ExponentialMeanReversion optimizer throws RootNotBracketedError when the
 * d*-equation has no root in the bracket.
 *
 */
TEST_F(TradingLevelsExponentialTest, entryLevelLowerExponentialNoRootTest) {
  // Declare and initialize test parameters.
  const double d_star = 1.236887;
  const double b_star = 1.4093;

  // Assert that the no-root bracket is reported.
  EXPECT_THROW(
      trading_levels.optimalEntryLower(d_star, b_star, rate, cost),
      RootNotBracketedError
  ) << "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntryLower did not "
       "throw RootNotBracketedError for a bracket with no root.";
}
/**
 * @test Solves every OU level on one parameter set, asserts the theory
 * ordering stop_loss < a* < d* < b*_L < b*, and asserts each solved value.
 *
 */
TEST_F(TradingLevelsTest, levelOrderingChainTest) {
  // Declare and initialize test parameters.
  const double stop_loss = 0.04;

  // Solve the full chain of levels.
  const double b_star = trading_levels.optimalExit(rate, cost);
  const double b_star_stop_loss =
      trading_levels.optimalExit(stop_loss, rate, cost);
  const double d_star =
      trading_levels.optimalEntry(b_star_stop_loss, stop_loss, rate, cost);
  const double a_star = trading_levels.optimalEntryLower(
      d_star, b_star_stop_loss, stop_loss, rate, cost
  );

  // Assert the ordering of the levels.
  EXPECT_LT(stop_loss, a_star)
      << "Lower entry level a* is not above the stop loss level.";
  EXPECT_LT(a_star, d_star)
      << "Entry level d* is not above the lower entry level a*.";
  EXPECT_LT(d_star, b_star_stop_loss)
      << "Stop-loss exit level is not above the entry level d*.";
  EXPECT_LT(b_star_stop_loss, b_star)
      << "Stop-loss exit level is not below the no-stop-loss exit level.";

  // Assert each solved level against the expected value.
  const double exit_tolerance = 5e-5;
  const double entry_tolerance = 2e-5;
  EXPECT_NEAR(
      b_star, expected::trading_levels_test::level_ordering_chain_test_b_star,
      exit_tolerance
  ) << "No-stop-loss exit level b* is not near the expected value.";
  EXPECT_NEAR(
      b_star_stop_loss,
      expected::trading_levels_test::level_ordering_chain_test_b_star_stop_loss,
      exit_tolerance
  ) << "Stop-loss exit level b*_L is not near the expected value.";
  EXPECT_NEAR(
      d_star, expected::trading_levels_test::level_ordering_chain_test_d_star,
      entry_tolerance
  ) << "Entry level d* is not near the expected value.";
  EXPECT_NEAR(
      a_star, expected::trading_levels_test::level_ordering_chain_test_a_star,
      entry_tolerance
  ) << "Lower entry level a* is not near the expected value.";
}
/**
 * @test Asserts that the exponential trading levels reject a stop-loss exit,
 * which has no mathematical definition for the exponential model.
 *
 */
TEST_F(
    TradingLevelsExponentialTest, exitLevelStopLossExponentialNotImplementedTest
) {
  // Declare and initialize test parameters.
  const double stop_loss = 0.04;

  // Assert that the method is not implemented.
  EXPECT_THROW(
      trading_levels.optimalExit(stop_loss, rate, cost), std::logic_error
  ) << "OrnsteinUhlenbeckTradingLevelsExponential::optimalExit with a stop "
       "loss did not throw.";
}
/**
 * @test Asserts that the exponential trading levels reject a stop-loss entry
 * level d*, which has no mathematical definition for the exponential model.
 *
 */
TEST_F(
    TradingLevelsExponentialTest,
    entryLevelStopLossExponentialNotImplementedTest
) {
  // Declare and initialize test parameters.
  const double stop_loss = 0.04;
  const double b_star = 1.5;

  // Assert that the method is not implemented.
  EXPECT_THROW(
      trading_levels.optimalEntry(b_star, stop_loss, rate, cost),
      std::logic_error
  ) << "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntry with a stop "
       "loss did not throw.";
}
/**
 * @test Asserts that the exponential trading levels reject a stop-loss lower
 * entry level a*, which has no mathematical definition for the exponential
 * model.
 *
 */
TEST_F(
    TradingLevelsExponentialTest,
    entryLevelLowerStopLossExponentialNotImplementedTest
) {
  // Declare and initialize test parameters.
  const double stop_loss = 0.04;
  const double b_star = 1.5;
  const double d_star = 1.0;

  // Assert that the method is not implemented.
  EXPECT_THROW(
      trading_levels.optimalEntryLower(d_star, b_star, stop_loss, rate, cost),
      std::logic_error
  ) << "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntryLower with a "
       "stop loss did not throw.";
}
/**
 * @test Solves every Ornstein-Uhlenbeck level on the fixture set and asserts
 * that each lies strictly inside the bracket its solver received.
 *
 */
TEST_F(TradingLevelsTest, solvedLevelsInsideBracketsTest) {
  const double stop_loss = 0.04;

  // Solve the levels in dependency order.
  const double b_star = trading_levels.optimalExit(rate, cost);
  const double b_star_stop_loss =
      trading_levels.optimalExit(stop_loss, rate, cost);
  const double d_star = trading_levels.optimalEntry(b_star, rate, cost);
  const double d_star_stop_loss =
      trading_levels.optimalEntry(b_star_stop_loss, stop_loss, rate, cost);
  const double a_star = trading_levels.optimalEntryLower(
      d_star_stop_loss, b_star_stop_loss, stop_loss, rate, cost
  );

  // Assert each level against the bracket its solver received.
  const double exit_lower = trading_levels.optimalExitLowerBound(rate, cost);
  const double exit_upper = trading_levels.optimalExitUpperBound();
  const double entry_lower = trading_levels.optimalEntryLowerBound();
  EXPECT_GT(b_star, exit_lower)
      << "Exit level b* is not above the exit lower bound.";
  EXPECT_LT(b_star, exit_upper)
      << "Exit level b* is not below the exit upper bound.";
  EXPECT_GT(b_star_stop_loss, exit_lower)
      << "Stop-loss exit level b*_L is not above the exit lower bound.";
  EXPECT_LT(b_star_stop_loss, exit_upper)
      << "Stop-loss exit level b*_L is not below the exit upper bound.";
  EXPECT_GT(d_star, entry_lower)
      << "Entry level d* is not above the entry lower bound.";
  EXPECT_LT(d_star, b_star) << "Entry level d* is not below the exit level b*.";
  EXPECT_GT(d_star_stop_loss, stop_loss)
      << "Stop-loss entry level d* is not above the stop loss.";
  EXPECT_LT(d_star_stop_loss, b_star_stop_loss)
      << "Stop-loss entry level d* is not below the stop-loss exit level b*_L.";
  EXPECT_GT(a_star, stop_loss)
      << "Lower entry level a* is not above the stop loss.";
  EXPECT_LT(a_star, d_star_stop_loss)
      << "Lower entry level a* is not below the stop-loss entry level d*.";
}
/**
 * @test Solves the 2 exponential levels on the fixture set and asserts that
 * each lies strictly inside the bracket its solver received.
 *
 */
TEST_F(TradingLevelsExponentialTest, solvedLevelsInsideBracketsTest) {
  // Solve the levels in dependency order.
  const double b_star = trading_levels.optimalExit(rate, cost);
  const double d_star = trading_levels.optimalEntry(b_star, rate, cost);

  // Assert each level against the bracket its solver received.
  EXPECT_GT(b_star, trading_levels.optimalExitLowerBound(rate, cost))
      << "Exponential exit level b* is not above the exit lower bound.";
  EXPECT_LT(b_star, trading_levels.optimalExitUpperBound())
      << "Exponential exit level b* is not below the exit upper bound.";
  EXPECT_GT(d_star, trading_levels.optimalEntryLowerBound())
      << "Exponential entry level d* is not above the entry lower bound.";
  EXPECT_LT(d_star, b_star)
      << "Exponential entry level d* is not below the exit level b*.";
}
