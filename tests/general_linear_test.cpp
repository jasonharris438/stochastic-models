#include "stochastic_models/sde/general_linear.h"

#include <gtest/gtest.h>
#include <stdexcept>
/**
 * @file
 * @brief Unit tests for the GeneralLinearModel class (mean/variance helpers).
 */

// Tests the return value of the getMean method.
TEST(GeneralLinearModelTest, GetMeanTest) {
  const float tolerance = 1e-5;
  const GeneralLinearModel model(-0.00143647, 10.4573);
  const double expected = 0.0;
  const double actual = model.getMean();
  EXPECT_EQ(expected, actual)
      << "GeneralLinearLikelihood getMean method returning invalid value.";
}

/**
 * @test Stationary variance sigma^2 / (2|mu|) exists only for mu < 0.
 */
TEST(GeneralLinearModelTest, GetUnconditionalVarianceTest) {
  const GeneralLinearModel model(-0.00143647, 10.4573);
  EXPECT_NEAR(model.getUnconditionalVariance(), 38063.83819015, 1e-3)
      << "getUnconditionalVariance not returning sigma^2 / (2|mu|).";
}

/**
 * @test No stationary variance exists for mu >= 0; the accessor must throw.
 */
TEST(GeneralLinearModelTest, GetUnconditionalVarianceThrowsForNonNegativeMu) {
  const GeneralLinearModel positive_mu(0.05, 0.4);
  EXPECT_THROW(positive_mu.getUnconditionalVariance(), std::domain_error)
      << "getUnconditionalVariance did not throw for mu > 0.";
  const GeneralLinearModel zero_mu(0.0, 0.4);
  EXPECT_THROW(zero_mu.getUnconditionalVariance(), std::domain_error)
      << "getUnconditionalVariance did not throw for mu = 0.";
}

/**
 * @test Conditional variance over a step t is
 * sigma^2 * (exp(2 mu t) - 1) / (2 mu).
 */
TEST(GeneralLinearModelTest, GetConditionalVarianceTest) {
  const GeneralLinearModel model(-0.00143647, 10.4573);
  EXPECT_NEAR(model.getConditionalVariance(1.0), 109.19818826, 1e-5)
      << "getConditionalVariance wrong at t = 1.";
  EXPECT_NEAR(model.getConditionalVariance(2.0), 218.08310690, 1e-5)
      << "getConditionalVariance wrong at t = 2.";
}

/**
 * @test Conditional variance at mu = 0 must take the sigma^2 * t limit
 * branch, and must agree with that limit near the mu = 1e-12 guard
 * threshold.
 */
TEST(GeneralLinearModelTest, GetConditionalVarianceZeroMuLimitTest) {
  const GeneralLinearModel zero_mu(0.0, 0.05);
  EXPECT_NEAR(zero_mu.getConditionalVariance(1.0), 0.0025, 1e-12)
      << "GeneralLinearModel conditional variance wrong at mu = 0, t = 1.";
  EXPECT_NEAR(zero_mu.getConditionalVariance(2.0), 0.005, 1e-12)
      << "GeneralLinearModel conditional variance wrong at mu = 0, t = 2.";

  // Just above the guard, so the closed-form branch is exercised at the
  // smallest mu it ever sees.
  const GeneralLinearModel near_zero_mu(1.1e-12, 0.05);
  EXPECT_NEAR(near_zero_mu.getConditionalVariance(1.0), 0.0025, 1e-12)
      << "GeneralLinearModel conditional variance discontinuous across the "
         "mu = 1e-12 guard threshold.";
}
/**
 * @test size = 0 must be rejected for a contract consistent with the OU model.
 */
TEST(GeneralLinearValidationTest, simulateRejectsZeroSize) {
  const GeneralLinearModel model(0.5, 0.1);
  EXPECT_THROW(model.Simulate(0.0, 0, 1), std::invalid_argument)
      << "Simulate accepted size == 0.";
}
