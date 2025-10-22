#include "stochastic_models/likelihood/general_linear_likelihood.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <cstdlib>
#include <gtest/gtest.h>
#include <unordered_map>

/**
 * @file
 * @brief Unit tests validating closed-form General Linear MLE calculations.
 *
 * These tests verify the values produced by `GeneralLinearLikelihood` match
 * precomputed reference results for a small sample series.
 */
TEST(GeneralLinearLikelihoodCalculateTest, ParameterTest) {
  const float tolerance = 1e-4;
  // Generate mock data.
  const std::vector<double> test_vec{1094.1, 1104.1, 1107.7, 1123.6, 1115.6,
                                     1112.7, 1118.4, 1116.9, 1127.9, 1153.2,
                                     1159.6, 1153.6, 1138.3, 1124.6, 1122.6,
                                     1134.,  1132.5, 1139.8, 1133.6, 1124.5};
  // Generate likelihood calculator and generate estimates.
  const GeneralLinearLikelihood likelihood;
  const GeneralLinearLikelihoodComponents components =
      likelihood.calculateComponents(test_vec);
  const GeneralLinearParameters params =
      likelihood.calculateParameters(components);

  // Expect equality for mu value.
  EXPECT_LE(abs(roundToDecimals(params.mu, 8) + 0.00143647), tolerance)
      << "GeneralLinearLikelihood not calculating correct value mu.";

  // Expect equality for sigma value.
  EXPECT_LE(abs(roundToDecimals(params.sigma, 8) - 10.4573), tolerance)
      << "GeneralLinearLikelihood not calculating correct value sigma.";
}
