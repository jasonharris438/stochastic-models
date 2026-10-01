#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h"
#include "support/expected_values.h"

#include <cmath>
#include <gtest/gtest.h>
#include <random>
#include <unordered_map>
#include <vector>

/**
 * @file
 * @brief Unit tests for Ornstein-Uhlenbeck likelihood component and parameter
 * calculations.
 */

/**
 * @test Tests that the Ornstein-Uhlenbeck likelihood class produces the
 * correct parameter estimates.
 *
 */
TEST(OrnsteinUhlenbeckLikelihoodCalculateTest, ParameterTest) {
  const double tolerance = 1e-9;
  // Generate mock data.
  const OrnsteinUhlenbeckLikelihoodComponents components = {4.0,   3.5,  4.125,
                                                            3.375, 3.25, 6};
  // Generate likelihood calculator and generate estimates.
  OrnsteinUhlenbeckLikelihood likelihood;
  const OrnsteinUhlenbeckParameters params =
      likelihood.calculateParameters(components);

  // Expect equality for mu value.
  EXPECT_NEAR(
      params.mu,
      expected::ornstein_uhlenbeck_likelihood_calculate_test::parameter_test_mu,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value mu.";
  // Expect equality for alpha value.
  EXPECT_NEAR(
      params.alpha,
      expected::ornstein_uhlenbeck_likelihood_calculate_test::
          parameter_test_alpha,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value alpha.";
  // Expect equality for sigma value.
  EXPECT_NEAR(
      params.sigma,
      expected::ornstein_uhlenbeck_likelihood_calculate_test::
          parameter_test_sigma,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value sigma.";
}

/**
 * @test Tests that the Ornstein-Uhlenbeck likelihood class returns the
 * correct updated parameter values.
 *
 */
TEST(OrnsteinUhlenbeckLikelihoodUpdateTest, ParameterTest) {
  const double tolerance = 1e-9;

  // Generate likelihood calculator and generate estimates.
  OrnsteinUhlenbeckLikelihood likelihood;
  const OrnsteinUhlenbeckLikelihoodComponents components = {4.0,   3.5,  4.125,
                                                            3.375, 3.25, 6};
  const OrnsteinUhlenbeckLikelihoodComponents updated_components =
      likelihood.updateComponents(components, 0.75, 1.0);

  const OrnsteinUhlenbeckParameters params =
      likelihood.calculateParameters(updated_components);

  // Expect equality for mu value.
  EXPECT_NEAR(
      params.mu,
      expected::ornstein_uhlenbeck_likelihood_update_test::parameter_test_mu,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value mu.";
  // Expect equality for alpha value.
  EXPECT_NEAR(
      params.alpha,
      expected::ornstein_uhlenbeck_likelihood_update_test::parameter_test_alpha,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value alpha.";
  // Expect equality for sigma value.
  EXPECT_NEAR(
      params.sigma,
      expected::ornstein_uhlenbeck_likelihood_update_test::parameter_test_sigma,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value sigma.";
}

/**
 * @test Tests that the Ornstein-Uhlenbeck likelihood class produces the
 * correct sufficient statistics.
 *
 */
TEST(OrnsteinUhlenbeckLikelihoodCalculateTest, ComponentsTest) {
  const double tolerance = 1e-12;
  // Generate mock data.
  const std::vector<double> test_vec{0.5, 0.25, 0.5, 0.75, 1.5, 1.0};
  // Generate likelihood calculator and generate estimates.
  OrnsteinUhlenbeckLikelihood likelihood;
  const OrnsteinUhlenbeckLikelihoodComponents components =
      likelihood.calculateComponents(test_vec);

  // Expect equality for the lead sum component.
  EXPECT_NEAR(
      components.lead_sum,
      expected::ornstein_uhlenbeck_likelihood_calculate_test::
          components_test_lead_sum,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
       "lead sum component.";
  // Expect equality for the lag sum component.
  EXPECT_NEAR(
      components.lag_sum,
      expected::ornstein_uhlenbeck_likelihood_calculate_test::
          components_test_lag_sum,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
       "lag sum component.";
  // Expect equality for the lead sum squared component.
  EXPECT_NEAR(
      components.lead_sum_squared,
      expected::ornstein_uhlenbeck_likelihood_calculate_test::
          components_test_lead_sum_squared,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lead "
       "sum squared component.";
  // Expect equality for the lag sum squared component.
  EXPECT_NEAR(
      components.lag_sum_squared,
      expected::ornstein_uhlenbeck_likelihood_calculate_test::
          components_test_lag_sum_squared,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lag "
       "sum squared component.";
  // Expect equality for the lead-lag sum product component.
  EXPECT_NEAR(
      components.lead_lag_sum_product,
      expected::ornstein_uhlenbeck_likelihood_calculate_test::
          components_test_lead_lag_sum_product,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
       "lead-lag sum product component.";
  // Expect equality for the number of observations component.
  EXPECT_EQ(
      components.n_obs, expected::ornstein_uhlenbeck_likelihood_calculate_test::
                            components_test_n_obs
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
       "number of observations component.";
}

/**
 * @test Tests that the Ornstein-Uhlenbeck likelihood class returns the
 * correct sufficient statistics after an update.
 *
 */
TEST(OrnsteinUhlenbeckLikelihoodUpdateTest, ComponentsTest) {
  const double tolerance = 1e-12;

  // Generate likelihood calculator and generate estimates.
  OrnsteinUhlenbeckLikelihood likelihood;
  const OrnsteinUhlenbeckLikelihoodComponents components = {4.0,   3.5,  4.125,
                                                            3.375, 3.25, 6};
  const OrnsteinUhlenbeckLikelihoodComponents updated_components =
      likelihood.updateComponents(components, 0.75, 1.0);

  // Expect equality for the lead sum component.
  EXPECT_NEAR(
      updated_components.lead_sum,
      expected::ornstein_uhlenbeck_likelihood_update_test::
          components_test_lead_sum,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lead "
       "sum component.";
  // Expect equality for the lag sum component.
  EXPECT_NEAR(
      updated_components.lag_sum,
      expected::ornstein_uhlenbeck_likelihood_update_test::
          components_test_lag_sum,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lag "
       "sum component.";
  // Expect equality for the lead sum squared component.
  EXPECT_NEAR(
      updated_components.lead_sum_squared,
      expected::ornstein_uhlenbeck_likelihood_update_test::
          components_test_lead_sum_squared,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lead "
       "sum squared component.";
  // Expect equality for the lag sum squared component.
  EXPECT_NEAR(
      updated_components.lag_sum_squared,
      expected::ornstein_uhlenbeck_likelihood_update_test::
          components_test_lag_sum_squared,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lag "
       "sum squared component.";
  // Expect equality for the lead-lag sum product component.
  EXPECT_NEAR(
      updated_components.lead_lag_sum_product,
      expected::ornstein_uhlenbeck_likelihood_update_test::
          components_test_lead_lag_sum_product,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
       "lead-lag sum product component.";
  // Expect equality for the number of observations component.
  EXPECT_EQ(
      updated_components.n_obs,
      expected::ornstein_uhlenbeck_likelihood_update_test::components_test_n_obs
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
       "number of observations component.";
}

// Data generated from the exact OU transition with known parameters must be
// recovered by the MLE. Tolerances are >= 6 standard errors at n = 20000.
TEST(OrnsteinUhlenbeckLikelihoodCalculateTest, RecoversKnownParameters) {
  const double mu = 1.3, alpha = 0.8, sigma = 0.5;
  const double a = std::exp(-alpha);
  const double transition_sd =
      sigma * std::sqrt((1 - std::exp(-2 * alpha)) / (2 * alpha));

  std::mt19937 gen(42);
  std::normal_distribution<double> norm(0.0, 1.0);
  std::vector<double> data{mu};
  for (int i = 0; i < 20000; i++) {
    data.push_back(mu + (data.back() - mu) * a + transition_sd * norm(gen));
  }

  OrnsteinUhlenbeckLikelihood likelihood;
  const OrnsteinUhlenbeckParameters params =
      likelihood.calculateParameters(likelihood.calculateComponents(data));

  EXPECT_NEAR(params.mu, mu, 0.03);
  EXPECT_NEAR(params.alpha, alpha, 0.10);
  EXPECT_NEAR(params.sigma, sigma, 0.02);
}

/**
 * @test An empty series must be rejected with InvalidNumberObservationsError
 * before any lead/lag iterator arithmetic runs.
 */
TEST(OuLikelihoodValidationTest, calculateComponentsRejectsEmptySeries) {
  const OrnsteinUhlenbeckLikelihood likelihood{};
  EXPECT_THROW(
      likelihood.calculateComponents({}), InvalidNumberObservationsError
  ) << "calculateComponents accepted an empty series.";
}

/**
 * @test A single-element series cannot form a lead/lag pair and must be
 * rejected with InvalidNumberObservationsError.
 */
TEST(OuLikelihoodValidationTest, calculateComponentsRejectsSingleObservation) {
  const OrnsteinUhlenbeckLikelihood likelihood{};
  EXPECT_THROW(
      likelihood.calculateComponents({1.0}), InvalidNumberObservationsError
  ) << "calculateComponents accepted a single-observation series.";
}

/**
 * @test calculateLeadSum must reject a series that cannot form a lead/lag
 * pair before any iterator arithmetic runs.
 */
TEST(OuLikelihoodValidationTest, calculateLeadSumRejectsShortSeries) {
  const OrnsteinUhlenbeckLikelihoodComponentCalculator calculator{};
  EXPECT_THROW(calculator.calculateLeadSum({}), InvalidNumberObservationsError)
      << "calculateLeadSum accepted an empty series.";
  EXPECT_THROW(
      calculator.calculateLeadSum({1.0}), InvalidNumberObservationsError
  ) << "calculateLeadSum accepted a single-observation series.";
}

/**
 * @test calculateLagSum previously advanced cend() by -1 on an empty series.
 * It must reject a short series.
 */
TEST(OuLikelihoodValidationTest, calculateLagSumRejectsShortSeries) {
  const OrnsteinUhlenbeckLikelihoodComponentCalculator calculator{};
  EXPECT_THROW(calculator.calculateLagSum({}), InvalidNumberObservationsError)
      << "calculateLagSum accepted an empty series.";
  EXPECT_THROW(
      calculator.calculateLagSum({1.0}), InvalidNumberObservationsError
  ) << "calculateLagSum accepted a single-observation series.";
}

/**
 * @test calculateLeadSumSquared must reject a short series.
 */
TEST(OuLikelihoodValidationTest, calculateLeadSumSquaredRejectsShortSeries) {
  const OrnsteinUhlenbeckLikelihoodComponentCalculator calculator{};
  EXPECT_THROW(
      calculator.calculateLeadSumSquared({}), InvalidNumberObservationsError
  ) << "calculateLeadSumSquared accepted an empty series.";
  EXPECT_THROW(
      calculator.calculateLeadSumSquared({1.0}), InvalidNumberObservationsError
  ) << "calculateLeadSumSquared accepted a single-observation series.";
}

/**
 * @test calculateLagSumSquared must reject a short series.
 */
TEST(OuLikelihoodValidationTest, calculateLagSumSquaredRejectsShortSeries) {
  const OrnsteinUhlenbeckLikelihoodComponentCalculator calculator{};
  EXPECT_THROW(
      calculator.calculateLagSumSquared({}), InvalidNumberObservationsError
  ) << "calculateLagSumSquared accepted an empty series.";
  EXPECT_THROW(
      calculator.calculateLagSumSquared({1.0}), InvalidNumberObservationsError
  ) << "calculateLagSumSquared accepted a single-observation series.";
}

/**
 * @test calculateLeadLagSumProduct must reject a short series.
 */
TEST(OuLikelihoodValidationTest, calculateLeadLagSumProductRejectsShortSeries) {
  const OrnsteinUhlenbeckLikelihoodComponentCalculator calculator{};
  EXPECT_THROW(
      calculator.calculateLeadLagSumProduct({}), InvalidNumberObservationsError
  ) << "calculateLeadLagSumProduct accepted an empty series.";
  EXPECT_THROW(
      calculator.calculateLeadLagSumProduct({1.0}),
      InvalidNumberObservationsError
  ) << "calculateLeadLagSumProduct accepted a single-observation series.";
}
