#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/numeric_utils/integration.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"
#include "stochastic_models/trading/exponential_mean_reversion.h"
#include "stochastic_models/trading/optimal_mean_reversion.h"
#include "stochastic_models/trading/trading_levels.h"
#include "stochastic_models/trading/trading_levels_exponential.h"
#include "stochastic_models/trading/trading_levels_params.h"

#include <gtest/gtest.h>
#include <stdexcept>
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
  const double stop_loss = 0.04;
  const double c = 0.02;
  const double r = 0.05;
  const double tolerance = 1e-4;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

  // Calculate b*.
  const double value = tradingLevels.optimalExit(stop_loss, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(value, 0.4551908, tolerance)
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
  const double tolerance = 5e-5;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

  // Calculate b*.
  const double value = tradingLevels.optimalExit(r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(value, 0.46683586991445225, tolerance)
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
  const double tolerance = 2e-4;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

  // Calculate b*.
  const double value = tradingLevels.optimalExit(r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(value, 1.40929954132349, tolerance)
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
  const double stop_loss = 0.04;
  const double d_star = 0.13093;
  const double b_star = 0.455191;
  const double c = 0.02;
  const double r = 0.05;
  const double tolerance = 1e-4;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

  // Calculate d*.
  const double value =
      tradingLevels.optimalEntryLower(d_star, b_star, stop_loss, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(value, 0.1076331, tolerance)
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
  const double stop_loss = 0.04;
  const double b_star = 0.455191;
  const double c = 0.02;
  const double r = 0.05;
  const double tolerance = 1e-4;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

  // Calculate d*.
  const double value = tradingLevels.optimalEntry(b_star, stop_loss, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(value, 0.1309298, tolerance)
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
  const double tolerance = 1e-4;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

  // Calculate d*.
  const double value = tradingLevels.optimalEntry(b_star, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(value, 0.1156811, tolerance)
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

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

  // Assert that the method is not implemented.
  ASSERT_THROW(tradingLevels.optimalEntry(b_star, r, c), std::invalid_argument)
      << "OrnsteinUhlenbeckTradingLevels::optimalEntry did not throw "
         "with invalid upper and lower bound values.";
}
/**
 * @test Asserts that the OU trading levels reject a no-stop-loss lower entry
 * level a*, which has no mathematical definition for the Ornstein-Uhlenbeck
 * model.
 *
 */
TEST(TradingLevelsTest, entryLevelLowerNotImplementedTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 8;
  const double mu = 0.3;
  const double sigma = 0.3;
  const double b_star = 0.4551908;
  const double d_star = 0.13093;
  const double c = 0.02;
  const double r = 0.05;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

  // Assert that the method is not implemented.
  EXPECT_THROW(
      tradingLevels.optimalEntryLower(d_star, b_star, r, c), std::logic_error
  ) << "OrnsteinUhlenbeckTradingLevels::optimalEntryLower without a stop "
       "loss did not throw.";
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
  const double tolerance = 1e-4;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

  // Calculate d*.
  const double value = tradingLevels.optimalEntry(b_star, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(value, 1.2368866, tolerance)
      << "Value produced by "
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
TEST(TradingLevelsTest, entryLevelLowerExponentialNoRootTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 5;
  const double mu = 1.3499;
  const double sigma = 0.15;
  const double d_star = 1.236887;
  const double b_star = 1.4093;
  const double c = 0.02;
  const double r = 0.05;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

  // Assert that the no-root bracket is reported.
  EXPECT_THROW(
      tradingLevels.optimalEntryLower(d_star, b_star, r, c),
      RootNotBracketedError
  ) << "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntryLower did not "
       "throw RootNotBracketedError for a bracket with no root.";
}
/**
 * @test Solves every OU level on one parameter set, asserts the theory
 * ordering stop_loss < a* < d* < b*_L < b*, and asserts each solved value.
 *
 */
TEST(TradingLevelsTest, levelOrderingChainTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 8;
  const double mu = 0.3;
  const double sigma = 0.3;
  const double stop_loss = 0.04;
  const double c = 0.02;
  const double r = 0.05;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

  // Solve the full chain of levels.
  const double b_star = tradingLevels.optimalExit(r, c);
  const double b_star_stop_loss = tradingLevels.optimalExit(stop_loss, r, c);
  const double d_star =
      tradingLevels.optimalEntry(b_star_stop_loss, stop_loss, r, c);
  const double a_star = tradingLevels.optimalEntryLower(
      d_star, b_star_stop_loss, stop_loss, r, c
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

  // Assert each solved level against the independently derived value.
  const double tolerance = 1e-4;
  EXPECT_NEAR(b_star, 0.4668359, tolerance)
      << "No-stop-loss exit level b* is not near the expected value.";
  EXPECT_NEAR(b_star_stop_loss, 0.4551908, tolerance)
      << "Stop-loss exit level b*_L is not near the expected value.";
  EXPECT_NEAR(d_star, 0.1309298, tolerance)
      << "Entry level d* is not near the expected value.";
  EXPECT_NEAR(a_star, 0.1076331, tolerance)
      << "Lower entry level a* is not near the expected value.";
}
/**
 * @test Asserts that the exponential trading levels reject a stop-loss exit,
 * which has no mathematical definition for the exponential model.
 *
 */
TEST(TradingLevelsTest, exitLevelStopLossExponentialNotImplementedTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 5;
  const double mu = 1.3499;
  const double sigma = 0.15;
  const double stop_loss = 0.04;
  const double c = 0.02;
  const double r = 0.05;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

  // Assert that the method is not implemented.
  EXPECT_THROW(tradingLevels.optimalExit(stop_loss, r, c), std::logic_error)
      << "OrnsteinUhlenbeckTradingLevelsExponential::optimalExit with a stop "
         "loss did not throw.";
}
/**
 * @test Asserts that the exponential trading levels reject a stop-loss entry
 * level d*, which has no mathematical definition for the exponential model.
 *
 */
TEST(TradingLevelsTest, entryLevelStopLossExponentialNotImplementedTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 5;
  const double mu = 1.3499;
  const double sigma = 0.15;
  const double stop_loss = 0.04;
  const double b_star = 1.5;
  const double c = 0.02;
  const double r = 0.05;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

  // Assert that the method is not implemented.
  EXPECT_THROW(
      tradingLevels.optimalEntry(b_star, stop_loss, r, c), std::logic_error
  ) << "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntry with a stop "
       "loss did not throw.";
}
/**
 * @test Asserts that the exponential trading levels reject a stop-loss lower
 * entry level a*, which has no mathematical definition for the exponential
 * model.
 *
 */
TEST(TradingLevelsTest, entryLevelLowerStopLossExponentialNotImplementedTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 5;
  const double mu = 1.3499;
  const double sigma = 0.15;
  const double stop_loss = 0.04;
  const double b_star = 1.5;
  const double d_star = 1.0;
  const double c = 0.02;
  const double r = 0.05;

  // Create trading levels instance to manage allocations.
  OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

  // Assert that the method is not implemented.
  EXPECT_THROW(
      tradingLevels.optimalEntryLower(d_star, b_star, stop_loss, r, c),
      std::logic_error
  ) << "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntryLower with a "
       "stop loss did not throw.";
}
