#include "stochastic_models/sde/ornstein_uhlenbeck.h"
#include "stochastic_models/trading/exponential_mean_reversion.h"
#include "stochastic_models/trading/optimal_mean_reversion.h"
#include "support/trading_fixtures.h"

#include <cstdlib>
#include <gtest/gtest.h>
#include <stdexcept>

namespace {

  class ExponentialMeanReversionTest
      : public test_support::TradingFixture<
            test_support::ornstein_uhlenbeck_set> {
  protected:
    HittingTimeOrnsteinUhlenbeck hitting_time_kernel{mu, alpha, sigma};
    ExponentialMeanReversion mean_reversion;
  };

} // namespace
/**
 * @test Tests the method ExponentialMeanReversion::b with a stop_loss argument
 * is not implemented
 *
 */
TEST_F(ExponentialMeanReversionTest, stopLossBNotImplementedTest) {
  // Declare and initialize model and test parameters.
  const double x = 0.2;
  const double stop_loss = 0.04;

  // Assert that the method is not implemented.
  ASSERT_THROW(
      mean_reversion.b(x, &hitting_time_kernel, stop_loss, rate, cost),
      std::logic_error
  ) << "OptimalMeanReversion::b with stop loss is not meant to be "
       "implemented.";
}
/**
 * @test Tests the method ExponentialMeanReversion::d with a stop_loss argument
 * is not implemented
 *
 */
TEST_F(ExponentialMeanReversionTest, stopLossDNotImplementedTest) {
  // Declare and initialize model and test parameters.
  const double x = 0.2;
  const double b_star = 0.466;
  const double stop_loss = 0.04;

  // Assert that the method is not implemented.
  ASSERT_THROW(
      mean_reversion.d(x, &hitting_time_kernel, b_star, stop_loss, rate, cost),
      std::logic_error
  ) << "OptimalMeanReversion::d with stop loss is not meant to be "
       "implemented.";
}
/**
 * @test Tests the method ExponentialMeanReversion::a with a stop_loss argument
 * is not implemented
 *
 */
TEST_F(ExponentialMeanReversionTest, stopLossANotImplementedTest) {
  // Declare and initialize model and test parameters.
  const double x = 0.2;
  const double b_star = 0.466;
  const double stop_loss = 0.04;

  // Assert that the method is not implemented.
  ASSERT_THROW(
      mean_reversion.a(x, &hitting_time_kernel, b_star, stop_loss, rate, cost),
      std::logic_error
  ) << "OptimalMeanReversion::a with stop loss is not meant to be "
       "implemented.";
}
/**
 * @test Tests the method ExponentialMeanReversion::V with a stop_loss argument
 * is not implemented
 *
 */
TEST_F(ExponentialMeanReversionTest, stopLossVNotImplementedTest) {
  // Declare and initialize model and test parameters.
  const double x = 0.2;
  const double b_star = 0.466;
  const double stop_loss = 0.04;

  // Assert that the method is not implemented.
  ASSERT_THROW(
      mean_reversion.V(&hitting_time_kernel, x, b_star, stop_loss, rate, cost),
      std::logic_error
  ) << "OptimalMeanReversion::V with stop loss is not meant to be "
       "implemented.";
}
/**
 * @test Tests the method ExponentialMeanReversion::instantaneousDifferential
 * with a stop_loss argument is not implemented
 *
 */
TEST_F(
    ExponentialMeanReversionTest,
    stopLossinstantaneousDifferentialNotImplementedTest
) {
  // Declare and initialize model and test parameters.
  const double x = 0.2;
  const double b_star = 0.466;
  const double stop_loss = 0.04;
  ModelFunc fn = funcIntegrateF;

  // Assert that the method is not implemented.
  ASSERT_THROW(
      mean_reversion.instantaneousDifferential(
          fn, &hitting_time_kernel, x, b_star, stop_loss, rate, cost
      ),
      std::logic_error
  ) << "OptimalMeanReversion::instantaneousDifferential with stop loss is "
       "not meant to be implemented.";
}
