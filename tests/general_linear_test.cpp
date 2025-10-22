#include "stochastic_models/numeric_utils/helpers.h"
#include "stochastic_models/sde/general_linear.h"

#include <gtest/gtest.h>
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

// Tests the return value of the getUnconditionalVariance method.
TEST(GeneralLinearModelTest, GetUnconditionalVarianceTest) {
  const float tolerance = 1e-5;
  const GeneralLinearModel model(-0.00143647, 10.4573);
  const double expected = 109.19818826;
  const double actual = model.getUnconditionalVariance();
  EXPECT_LE(abs(roundToDecimals(actual, 8) - expected), tolerance)
      << "GeneralLinearLikelihood getUnconditionalVariance method returning "
         "invalid value.";
}

// Tests the return value of the getConditionalVariance method.
TEST(GeneralLinearModelTest, GetConditionalVarianceTest) {
  const float tolerance = 1e-5;
  const GeneralLinearModel model(-0.00143647, 10.4573);
  const double expected = 20.95971157;
  const double actual = model.getConditionalVariance();
  EXPECT_LE(abs(roundToDecimals(actual, 8) - expected), tolerance)
      << "GeneralLinearLikelihood getConditionalVariance method returning "
         "invalid value.";
}
