#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/likelihood/general_linear_likelihood.h"
#include "support/expected_values.h"
#include "support/general_linear_fixtures.h"

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
  const double tolerance = 1e-9;
  // Generate mock data.
  const std::vector<double> test_vec = test_support::generalLinearSeries();
  // Generate likelihood calculator and generate estimates.
  const GeneralLinearLikelihood likelihood;
  const GeneralLinearLikelihoodComponents components =
      likelihood.calculateComponents(test_vec);
  const GeneralLinearParameters params =
      likelihood.calculateParameters(components);

  // Expect equality for mu value.
  EXPECT_NEAR(
      params.mu,
      expected::general_linear_likelihood_calculate_test::parameter_test_mu,
      tolerance
  ) << "GeneralLinearLikelihood not calculating correct value mu.";

  // Expect equality for sigma value.
  EXPECT_NEAR(
      params.sigma,
      expected::general_linear_likelihood_calculate_test::parameter_test_sigma,
      tolerance
  ) << "GeneralLinearLikelihood not calculating correct value sigma.";
}

/**
 * @test An empty series must be rejected before the lead/lag helpers run —
 * previously the size guard sat in calculateSquaredError, which is reached
 * only after the unguarded iterator arithmetic.
 */
TEST(
    GeneralLinearLikelihoodValidationTest, calculateComponentsRejectsEmptySeries
) {
  const GeneralLinearLikelihood likelihood{};
  EXPECT_THROW(
      likelihood.calculateComponents({}), InvalidNumberObservationsError
  ) << "calculateComponents accepted an empty series.";
}

/**
 * @test A single-element series cannot form a lead/lag pair and must be
 * rejected with InvalidNumberObservationsError.
 */
TEST(
    GeneralLinearLikelihoodValidationTest,
    calculateComponentsRejectsSingleObservation
) {
  const GeneralLinearLikelihood likelihood{};
  EXPECT_THROW(
      likelihood.calculateComponents({1.0}), InvalidNumberObservationsError
  ) << "calculateComponents accepted a single-observation series.";
}

/**
 * @test calculateLeadLagInnerProduct previously advanced cend() by -1 on an
 * empty series. It must reject a short series.
 */
TEST(
    GeneralLinearLikelihoodValidationTest,
    calculateLeadLagInnerProductRejectsShortSeries
) {
  const GeneralLinearLikelihoodComponentCalculator calculator{};
  EXPECT_THROW(
      calculator.calculateLeadLagInnerProduct({}),
      InvalidNumberObservationsError
  ) << "calculateLeadLagInnerProduct accepted an empty series.";
  EXPECT_THROW(
      calculator.calculateLeadLagInnerProduct({1.0}),
      InvalidNumberObservationsError
  ) << "calculateLeadLagInnerProduct accepted a single-observation series.";
}

/**
 * @test calculateLagSquared must reject a short series.
 */
TEST(
    GeneralLinearLikelihoodValidationTest, calculateLagSquaredRejectsShortSeries
) {
  const GeneralLinearLikelihoodComponentCalculator calculator{};
  EXPECT_THROW(
      calculator.calculateLagSquared({}), InvalidNumberObservationsError
  ) << "calculateLagSquared accepted an empty series.";
  EXPECT_THROW(
      calculator.calculateLagSquared({1.0}), InvalidNumberObservationsError
  ) << "calculateLagSquared accepted a single-observation series.";
}

/**
 * @test calculateSquaredError must reject an empty series and a one-element
 * series with the library's observation-count exception.
 */
TEST(
    GeneralLinearLikelihoodValidationTest,
    calculateSquaredErrorRejectsShortSeries
) {
  const GeneralLinearLikelihoodComponentCalculator calculator{};
  EXPECT_THROW(
      calculator.calculateSquaredError({}, 0.0), InvalidNumberObservationsError
  ) << "calculateSquaredError accepted an empty series.";
  EXPECT_THROW(
      calculator.calculateSquaredError({1.0}, 0.0),
      InvalidNumberObservationsError
  ) << "calculateSquaredError accepted a one-element series.";
}

/**
 * @test calculateLagSquared sums the squares of the lagged values, so the
 * last observation is excluded.
 */
TEST(GeneralLinearLikelihoodCalculateTest, LagSquaredExcludesLastObservation) {
  const GeneralLinearLikelihoodComponentCalculator calculator{};
  EXPECT_DOUBLE_EQ(calculator.calculateLagSquared({1.0, 2.0, 3.0}), 5.0)
      << "calculateLagSquared did not sum the lagged squares 1 + 4.";
}
