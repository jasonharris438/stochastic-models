#include "stochastic_models/likelihood/general_linear_likelihood.h"
#include "stochastic_models/likelihood/general_linear_online.h"
#include "support/expected_values.h"
#include "support/general_linear_fixtures.h"

#include <cmath>
#include <gtest/gtest.h>
#include <vector>

/**
 * @file
 * @brief Unit tests for the online General Linear MLE updater.
 *
 * These tests exercise single-step updates for `mu` and `sigma` using a
 * fixed time series and assert values computed by
 * `tools/derivations/expected_values.py`.
 */

/**
 * @test Tests that the GeneralLinearUpdater.updateMu method returns the
 * correct value.
 *
 */
TEST(GeneralLinearOnlineTest, UpdateMuTest) {
  // Create test data.
  const std::vector<double> test_vec = test_support::generalLinearSeries();
  const double tolerance = 1e-9;
  const double new_observation = test_support::general_linear_next_observation;
  const double initial_observation = 1124.5;

  // Calculate expected value tracking components.
  const GeneralLinearLikelihood likelihood;
  const GeneralLinearLikelihoodComponents components =
      likelihood.calculateComponents(test_vec);
  const GeneralLinearParameters params =
      likelihood.calculateParameters(components);

  GeneralLinearUpdater updater{components, params};

  // Calculate the new parameters.
  const GeneralLinearParameters actual =
      updater.updateState(new_observation, initial_observation);

  const double expected = expected::general_linear_online_test::update_mu_test;
  EXPECT_NEAR(actual.mu, expected, tolerance)
      << "GeneralLinearUpdater updateMu method returning invalid value.";
}

/**
 * @test Tests that the GeneralLinearUpdater.updateSigma method returns the
 * correct value.
 *
 */
TEST(GeneralLinearOnlineTest, UpdateSigmaTest) {
  // Create test data.
  const std::vector<double> test_vec = test_support::generalLinearSeries();
  const double tolerance = 1e-9;
  const double new_observation = test_support::general_linear_next_observation;
  const double initial_observation = 1124.5;

  // Calculate expected value tracking components.
  const GeneralLinearLikelihood likelihood;
  const GeneralLinearLikelihoodComponents components =
      likelihood.calculateComponents(test_vec);
  const GeneralLinearParameters params =
      likelihood.calculateParameters(components);

  GeneralLinearUpdater updater{components, params};

  // Calculate the new parameters.
  const GeneralLinearParameters actual =
      updater.updateState(new_observation, initial_observation);

  const double expected =
      expected::general_linear_online_test::update_sigma_test;
  EXPECT_NEAR(actual.sigma, expected, tolerance)
      << "GeneralLinearUpdater updateSigma method returning invalid value.";
}

/**
 * @test 1 online update from the 20-point series equals a batch fit of the
 * 21-point series, so the squared-error update is exact.
 */
TEST(GeneralLinearOnlineTest, OneUpdateEqualsBatchRecalculation) {
  const double parameter_tolerance = 1e-9;
  const double relative_component_tolerance = 1e-12;
  const GeneralLinearLikelihood likelihood;
  const std::vector<double> series = test_support::generalLinearSeries();
  std::vector<double> extended_series = series;
  extended_series.push_back(test_support::general_linear_next_observation);

  const GeneralLinearLikelihoodComponents batch_components =
      likelihood.calculateComponents(extended_series);
  const GeneralLinearParameters batch_parameters =
      likelihood.calculateParameters(batch_components);

  const GeneralLinearLikelihoodComponents online_components =
      likelihood.updateComponents(
          likelihood.calculateComponents(series),
          test_support::general_linear_next_observation, series.back()
      );
  const GeneralLinearParameters online_parameters =
      likelihood.calculateParameters(online_components);

  EXPECT_NEAR(
      online_components.lag_squared, batch_components.lag_squared,
      relative_component_tolerance * std::abs(batch_components.lag_squared)
  ) << "Online lag_squared differs from the batch value.";
  EXPECT_NEAR(
      online_components.lead_lag_inner_product,
      batch_components.lead_lag_inner_product,
      relative_component_tolerance *
          std::abs(batch_components.lead_lag_inner_product)
  ) << "Online lead_lag_inner_product differs from the batch value.";
  EXPECT_NEAR(
      online_components.squared_error, batch_components.squared_error,
      relative_component_tolerance * std::abs(batch_components.squared_error)
  ) << "Online squared_error differs from the batch value.";
  EXPECT_EQ(online_components.n_obs, batch_components.n_obs)
      << "Online n_obs differs from the batch value.";
  EXPECT_NEAR(online_parameters.mu, batch_parameters.mu, parameter_tolerance)
      << "Online mu differs from the batch value.";
  EXPECT_NEAR(
      online_parameters.sigma, batch_parameters.sigma, parameter_tolerance
  ) << "Online sigma differs from the batch value.";
}
