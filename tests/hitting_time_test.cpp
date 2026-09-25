#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/hitting_times/hitting_time_density.h"
#include "stochastic_models/hitting_times/hitting_time_ornstein_uhlenbeck.h"
#include "stochastic_models/numeric_utils/integration.h"
#include "support/expected_values.h"

#include <gtest/gtest.h>
/**
 * @test Tests the output of the hittingTimeDensity function and asserts that
 * it is near the expected value.
 *
 */
TEST(HittingTimeDensityTest, hittingTimeDensityOutputTest) {
  // Declare and initialize model and test parameters.
  double alpha = 0.0045;
  double mu = 0.998;
  double sigma = 0.0038;
  double first = 1.04;
  double second = 1;
  double x = 1.02;
  const double tolerance = 1e-6;
  // Create core model instance and declare function to use.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel =
      HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  ModelFunc fn = &integrateHittingTimeDensity;

  // Calculate the hitting time density.
  const double value =
      hittingTimeDensity(x, fn, &hitting_time_kernel, first, second);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value,
      expected::hitting_time_density_test::hitting_time_density_output_test,
      tolerance
  ) << "The value of the hitting time density is not equal to the expected "
       "value.";
}

/**
 * @test Tests that the hittingTimeDensity function throws the correct exception
 * when an incorrectly specified void pointer is provided.
 *
 */
TEST(HittingTimeDensityTest, CatchErrorTest) {
  double first = 1.04;
  double second = 1;
  double x = 1.02;
  std::string hitting_time_kernel = "This is not a model";
  ModelFunc fn = &integrateHittingTimeDensity;
  EXPECT_THROW(
      hittingTimeDensity(x, fn, &hitting_time_kernel, first, second),
      std::runtime_error
  ) << "hittingTimeDensity function must throw a runtime_error when a void "
       "pointer to the wrong type is provided.";
}

// Integrand that makes both density integrals exactly zero.
double zeroIntegrand(double, void*) {
  return 0.0;
}

/**
 * @test Tests that hittingTimeDensity throws ZeroDivError when the
 * denominator integral is zero.
 *
 */
TEST(HittingTimeDensityValidationTest, ZeroDenominatorThrowTest) {
  double first = 1.04;
  double second = 1;
  double x = 1.02;

  ModelFunc fn = &zeroIntegrand;

  EXPECT_THROW(hittingTimeDensity(x, fn, nullptr, first, second), ZeroDivError)
      << "hittingTimeDensity did not throw ZeroDivError for a zero "
         "denominator integral.";
}
