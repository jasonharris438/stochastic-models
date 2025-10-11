#include "stochastic_models/numeric_utils/helpers.h"
#include "stochastic_models/numeric_utils/integration.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"
#include "stochastic_models/trading/optimal_mean_reversion.h"

#include <cstdlib>
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
  const float tolerance = 1e-5;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate F(x;r).
  const double value = mean_reversion.F(hitting_time_kernel, x, r, c);
  delete hitting_time_kernel;

  // Assert that the value is near the expected value.
  EXPECT_LE(abs(roundToDecimals(value, 8) - 5.80482400), tolerance)
      << "Value produced by OptimalMeanReversion::F is not equal to the "
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
  const float tolerance = 2e-4;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate G(x;r).
  const double value = mean_reversion.G(hitting_time_kernel, x, r, c);
  delete hitting_time_kernel;

  // Assert that the value is near the expected value.
  EXPECT_LE(abs(roundToDecimals(value, 8) - 163.617), tolerance)
      << "Value produced by OptimalMeanReversion::G is not equal to the "
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
  const float tolerance = 1e-3;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate b().
  const double value = mean_reversion.b(guess, hitting_time_kernel, r, c);
  delete hitting_time_kernel;

  // Assert that the value is near the expected value.
  EXPECT_LE(abs(roundToDecimals(value, 8) - 135.511), tolerance)
      << "Value produced by OptimalMeanReversion::b is not equal to the "
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
  const float tolerance = 1e-4;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate b().
  const double value =
      mean_reversion.b(guess, hitting_time_kernel, stop_loss, r, c);
  delete hitting_time_kernel;

  // Assert that the value is near the expected value.
  EXPECT_LE(abs(roundToDecimals(value, 2) - 187.81), tolerance)
      << "Value produced by OptimalMeanReversion::b with stop loss level is "
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
  const float tolerance = 1e-3;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate d().
  const double value =
      mean_reversion.d(guess, hitting_time_kernel, b_star, r, c);
  delete hitting_time_kernel;

  // Assert that the value is near the expected value.
  EXPECT_LE(abs(roundToDecimals(value, 8) + 129.094), tolerance)
      << "Value produced by OptimalMeanReversion::d is not equal to the "
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
  const double stop_loss = -0.3;
  const double b_star = 0.466836;
  const double c = 0.02;
  const double r = 0.05;
  const float tolerance = 1e-3;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate d().
  const double value =
      mean_reversion.d(guess, hitting_time_kernel, b_star, stop_loss, r, c);
  delete hitting_time_kernel;

  // Assert that the value is near the expected value.
  EXPECT_LE(abs(roundToDecimals(value, 8) + 129.094), tolerance)
      << "Value produced by OptimalMeanReversion::d is not equal to the "
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
  const float tolerance = 1e-3;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate a().
  const double value =
      mean_reversion.a(guess, hitting_time_kernel, b_star, stop_loss, r, c);
  delete hitting_time_kernel;

  // Assert that the value is near the expected value.
  EXPECT_LE(abs(roundToDecimals(value, 8) + 132.55), tolerance)
      << "Value produced by OptimalMeanReversion::a is not equal to the "
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
  const float tolerance = 1e-5;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate V(x).
  const double value = mean_reversion.V(hitting_time_kernel, x, b_star, r, c);
  delete hitting_time_kernel;

  // Assert that the value is near the expected value.
  EXPECT_LE(abs(roundToDecimals(value, 8) - 0.53), tolerance)
      << "Value produced by OptimalMeanReversion::V is not equal to the "
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
  const float tolerance = 1e-5;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate V(x).
  const double value = mean_reversion.V(hitting_time_kernel, x, b_star, r, c);
  delete hitting_time_kernel;

  // Assert that the value is near the expected value.
  EXPECT_LE(abs(roundToDecimals(value, 8) - 0.403389), tolerance)
      << "Value produced by OptimalMeanReversion::V is not equal to the "
         "expected value when x is below b*.";
}
/**
 * @test Tests the output of the OptimalMeanReversion::V method when a stop loss
 * is provided and the x value is above b*. Asserts that it is near the
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
  const float tolerance = 1e-5;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate V(x).
  const double value =
      mean_reversion.V(hitting_time_kernel, x, b_star, stop_loss, r, c);
  delete hitting_time_kernel;

  // Assert that the value is near the expected value.
  EXPECT_LE(abs(roundToDecimals(value, 8) - 0.55), tolerance)
      << "Value produced by OptimalMeanReversion::V is not equal to the "
         "expected value when x is above b* and a stop loss is provided.";
}
/**
 * @test Tests the output of the OptimalMeanReversion::V method when a stop loss
 * is provided and the x value is below b*. Asserts that it is near the
 * expected value.
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
  const double x = 0.15;
  const float tolerance = 1e-5;

  // Create core model and optimal mean reversion instances.
  HittingTimeOrnsteinUhlenbeck* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  OptimalMeanReversion mean_reversion;

  // Calculate V(x).
  const double value = mean_reversion.V(hitting_time_kernel, x, b_star, r, c);
  delete hitting_time_kernel;

  // Assert that the value is near the expected value.
  EXPECT_LE(abs(roundToDecimals(value, 8) - 0.508892), tolerance)
      << "Value produced by OptimalMeanReversion::V is not equal to the "
         "expected value when x is below b* and a stop loss is provided.";
}
