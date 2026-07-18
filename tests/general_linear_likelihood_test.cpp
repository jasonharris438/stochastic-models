#include "stochastic_models/likelihood/general_linear_likelihood.h"
#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/numeric_utils/helpers.h"

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
  EXPECT_NEAR(roundToDecimals(params.mu, 8), -0.00143647, tolerance)
      << "GeneralLinearLikelihood not calculating correct value mu.";

  // Expect equality for sigma value.
  EXPECT_NEAR(roundToDecimals(params.sigma, 8), 10.4573, tolerance)
      << "GeneralLinearLikelihood not calculating correct value sigma.";
}

/**
 * @test An empty series must be rejected before the lead/lag helpers run —
 * previously the size guard sat in calculateSquaredError, which is reached
 * only after the unguarded iterator arithmetic.
 */
TEST(GeneralLinearLikelihoodValidationTest, calculateComponentsRejectsEmptySeries) {
  const GeneralLinearLikelihood likelihood{};
  EXPECT_THROW(
      likelihood.calculateComponents({}), InvalidNumberObservationsError
  ) << "calculateComponents accepted an empty series.";
}

/**
 * @test A single-element series cannot form a lead/lag pair and must be
 * rejected with InvalidNumberObservationsError.
 */
TEST(GeneralLinearLikelihoodValidationTest, calculateComponentsRejectsSingleObservation) {
  const GeneralLinearLikelihood likelihood{};
  EXPECT_THROW(
      likelihood.calculateComponents({1.0}), InvalidNumberObservationsError
  ) << "calculateComponents accepted a single-observation series.";
}
