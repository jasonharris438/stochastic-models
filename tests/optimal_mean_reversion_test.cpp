#include "stochastic_models/numeric_utils/integration.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"
#include "stochastic_models/trading/optimal_mean_reversion.h"
#include "support/expected_values.h"

#include <gtest/gtest.h>
/**
 * @test Tests the output of the OptimalMeanReversion::F method and asserts that
 * it is near the expected value.
 *
 */
TEST(OptimalMeanReversionTest, methodFOutputTest) {
  // Declare and initialize model and test parameters.
  const double sigma = 0.003;
  const double mu = 0.995;
  const double alpha = 0.02;
  const double r = 0.05;
  const double c = 0.001;
  const double x = 1.01;
  const double tolerance = 1e-5;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate F(x;r).
  const double value = mean_reversion.F(&hitting_time_kernel, x, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value, expected::optimal_mean_reversion_test::method_f_output_test,
      tolerance
  ) << "Value produced by OptimalMeanReversion::F is not equal to the "
       "expected value.";
}
/**
 * @test Tests the output of the OptimalMeanReversion::G method and asserts that
 * it is near the expected value.
 *
 */
TEST(OptimalMeanReversionTest, methodGOutputTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 8;
  const double mu = 0.3;
  const double sigma = 0.3;
  const double x = 0.2;
  const double c = 0.02;
  const double r = 0.05;
  const double tolerance = 2e-4;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate G(x;r).
  const double value = mean_reversion.G(&hitting_time_kernel, x, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value, expected::optimal_mean_reversion_test::method_g_output_test,
      tolerance
  ) << "Value produced by OptimalMeanReversion::G is not equal to the "
       "expected value.";
}
/**
 * @test Tests the output of the OptimalMeanReversion::b method and asserts that
 * it is near the expected value.
 *
 */
TEST(OptimalMeanReversionTest, methodBOutputTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 8;
  const double mu = 0.3;
  const double sigma = 0.3;
  const double guess = 0.4;
  const double c = 0.02;
  const double r = 0.05;
  const double tolerance = 1e-3;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate b().
  const double value = mean_reversion.b(guess, &hitting_time_kernel, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value, expected::optimal_mean_reversion_test::method_b_output_test,
      tolerance
  ) << "Value produced by OptimalMeanReversion::b is not equal to the "
       "expected value.";
}
/**
 * @test Tests the output of the OptimalMeanReversion::b method when a stop loss
 * level is provided and asserts that it is near the expected value.
 *
 */
TEST(OptimalMeanReversionTest, methodBStopLossOutputTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 8;
  const double mu = 0.3;
  const double stop_loss = 0.2;
  const double sigma = 0.3;
  const double guess = 0.28;
  const double c = 0.02;
  const double r = 0.05;
  const double tolerance = 5e-3;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate b().
  const double value =
      mean_reversion.b(guess, &hitting_time_kernel, stop_loss, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::optimal_mean_reversion_test::method_b_stop_loss_output_test,
      tolerance
  ) << "Value produced by OptimalMeanReversion::b with stop loss level is "
       "not equal to the expected value.";
}
/**
 * @test Tests the output of the OptimalMeanReversion::d method and asserts that
 * it is near the expected value.
 *
 */
TEST(OptimalMeanReversionTest, methodDOutputTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 8;
  const double mu = 0.3;
  const double sigma = 0.3;
  const double guess = 0.4;
  const double b_star = 0.46683583;
  const double c = 0.02;
  const double r = 0.05;
  const double tolerance = 1e-3;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate d().
  const double value =
      mean_reversion.d(guess, &hitting_time_kernel, b_star, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value, expected::optimal_mean_reversion_test::method_d_output_test,
      tolerance
  ) << "Value produced by OptimalMeanReversion::d is not equal to the "
       "expected value.";
}
/**
 * @test Tests the output of the OptimalMeanReversion::d method when a stop loss
 * is provided and asserts that it is near the expected value.
 *
 */
TEST(OptimalMeanReversionTest, methodDStopLossOutputTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 8;
  const double mu = 0.3;
  const double sigma = 0.3;
  const double guess = 0.4;
  const double stop_loss = 0.1;
  const double b_star = 0.466836;
  const double c = 0.02;
  const double r = 0.05;
  const double tolerance = 1e-3;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate d().
  const double value =
      mean_reversion.d(guess, &hitting_time_kernel, b_star, stop_loss, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::optimal_mean_reversion_test::method_d_stop_loss_output_test,
      tolerance
  ) << "Value produced by OptimalMeanReversion::d is not equal to the "
       "expected value when a stop loss is provided.";
}
/**
 * @test Tests the output of the OptimalMeanReversion::a method when a stop loss
 * is provided and asserts that it is near the expected value.
 *
 */
TEST(OptimalMeanReversionTest, methodAStopLossOutputTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 8;
  const double mu = 0.3;
  const double sigma = 0.3;
  const double guess = 0.4;
  const double stop_loss = -0.3;
  const double b_star = 0.466836;
  const double c = 0.02;
  const double r = 0.05;
  const double tolerance = 1e-3;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate a().
  const double value =
      mean_reversion.a(guess, &hitting_time_kernel, b_star, stop_loss, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::optimal_mean_reversion_test::method_a_stop_loss_output_test,
      tolerance
  ) << "Value produced by OptimalMeanReversion::a is not equal to the "
       "expected value when a stop loss is provided.";
}
/**
 * @test Tests the output of the OptimalMeanReversion::V method when the x value
 * is above b_star. Asserts that it is near the expected value.
 *
 */
TEST(OptimalMeanReversionTest, methodAboveVOutputTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 8;
  const double mu = 0.3;
  const double sigma = 0.3;
  const double b_star = 0.466836;
  const double c = 0.02;
  const double r = 0.05;
  const double x = 0.55;
  const double tolerance = 1e-9;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate V(x).
  const double value = mean_reversion.V(&hitting_time_kernel, x, b_star, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value, expected::optimal_mean_reversion_test::method_above_v_output_test,
      tolerance
  ) << "Value produced by OptimalMeanReversion::V is not equal to the "
       "expected value when x is above b*.";
}
/**
 * @test Tests the output of the OptimalMeanReversion::V method when the x value
 * is below b_star. Asserts that it is near the expected value.
 *
 */
TEST(OptimalMeanReversionTest, methodBelowVOutputTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 8;
  const double mu = 0.3;
  const double sigma = 0.3;
  const double b_star = 0.466836;
  const double c = 0.02;
  const double r = 0.05;
  const double x = 0.15;
  const double tolerance = 1e-6;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate V(x).
  const double value = mean_reversion.V(&hitting_time_kernel, x, b_star, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value, expected::optimal_mean_reversion_test::method_below_v_output_test,
      tolerance
  ) << "Value produced by OptimalMeanReversion::V is not equal to the "
       "expected value when x is below b*.";
}
/**
 * @test Tests that the OptimalMeanReversion::V method with a stop loss
 * returns the payoff x - c when x is above b*. Asserts that it is near the
 * expected value.
 *
 */
TEST(OptimalMeanReversionTest, methodAboveVStopLossOutputTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 16.6677;
  const double mu = 0.5388;
  const double sigma = 0.1599;
  const double stop_loss = 0.4834;
  const double b_star = 0.567304;
  const double c = 0.05;
  const double r = 0.05;
  const double x = 0.6;
  const double tolerance = 1e-9;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate V(x).
  const double value =
      mean_reversion.V(&hitting_time_kernel, x, b_star, stop_loss, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::optimal_mean_reversion_test::
          method_above_v_stop_loss_output_test,
      tolerance
  ) << "Value produced by OptimalMeanReversion::V is not equal to the "
       "expected value when x is above b* and a stop loss is provided.";
}
/**
 * @test Tests the output of the OptimalMeanReversion::V method when a stop loss
 * is provided and the x value is between the stop loss and b*. Asserts that it
 * is near the expected value.
 *
 */
TEST(OptimalMeanReversionTest, methodBelowVStopLossOutputTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 16.6677;
  const double mu = 0.5388;
  const double sigma = 0.1599;
  const double stop_loss = 0.4834;
  const double b_star = 0.567304;
  const double c = 0.05;
  const double r = 0.05;
  const double x = 0.5;
  const double tolerance = 1e-6;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate V(x).
  const double value =
      mean_reversion.V(&hitting_time_kernel, x, b_star, stop_loss, r, c);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::optimal_mean_reversion_test::
          method_below_v_stop_loss_output_test,
      tolerance
  ) << "Value produced by OptimalMeanReversion::V is not equal to the "
       "expected value when x is between the stop loss and b*.";
}

/**
 * @test Tests that the stop-loss value function returns the expected value
 * 1e-6 inside each boundary of its continuation region. There it lies within
 * 4e-6 of the boundary payoff.
 *
 */
TEST(OptimalMeanReversionTest, methodVStopLossBoundaryTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 16.6677;
  const double mu = 0.5388;
  const double sigma = 0.1599;
  const double stop_loss = 0.4834;
  const double b_star = 0.567304;
  const double c = 0.05;
  const double r = 0.05;
  const double offset = 1e-6;
  const double tolerance = 1e-6;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Evaluate V(x) just inside each boundary of the continuation region.
  const double at_exit = mean_reversion.V(
      &hitting_time_kernel, b_star - offset, b_star, stop_loss, r, c
  );
  const double at_stop_loss = mean_reversion.V(
      &hitting_time_kernel, stop_loss + offset, b_star, stop_loss, r, c
  );

  // Assert that both values are near the expected values.
  EXPECT_NEAR(
      at_exit,
      expected::optimal_mean_reversion_test::
          method_v_stop_loss_boundary_test_at_exit,
      tolerance
  ) << "Value produced by OptimalMeanReversion::V just below b* with a stop "
       "loss is not the expected value.";
  EXPECT_NEAR(
      at_stop_loss,
      expected::optimal_mean_reversion_test::
          method_v_stop_loss_boundary_test_at_stop_loss,
      tolerance
  ) << "Value produced by OptimalMeanReversion::V just above the stop loss "
       "is not the expected value.";
}
