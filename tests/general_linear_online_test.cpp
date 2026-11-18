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
 * These tests exercise online updates for `mu` and `sigma` using a fixed
 * time series and assert values computed by
 * `tools/derivations/expected_values.py`.
 */

namespace {

  /**
   * @brief Expects online components and their parameters to equal a batch
   * fit of `series`: sums within 1e-12 relative, the count exactly, and
   * parameters within 1e-9. `stream_description` begins each failure message.
   */
  void expectOnlineEqualsBatch(
      const GeneralLinearLikelihood& likelihood,
      const GeneralLinearLikelihoodComponents& online_components,
      const std::vector<double>& series,
      const char* stream_description
  ) {
    const double parameter_tolerance = 1e-9;
    const double relative_component_tolerance = 1e-12;
    const GeneralLinearLikelihoodComponents batch_components =
        likelihood.calculateComponents(series);
    const GeneralLinearParameters batch_parameters =
        likelihood.calculateParameters(batch_components);
    const GeneralLinearParameters online_parameters =
        likelihood.calculateParameters(online_components);

    EXPECT_NEAR(
        online_components.lag_squared, batch_components.lag_squared,
        relative_component_tolerance * std::abs(batch_components.lag_squared)
    ) << stream_description
      << " differs from the batch fit in lag_squared.";
    EXPECT_NEAR(
        online_components.lead_lag_inner_product,
        batch_components.lead_lag_inner_product,
        relative_component_tolerance *
            std::abs(batch_components.lead_lag_inner_product)
    ) << stream_description
      << " differs from the batch fit in lead_lag_inner_product.";
    EXPECT_NEAR(
        online_components.squared_error, batch_components.squared_error,
        relative_component_tolerance * std::abs(batch_components.squared_error)
    ) << stream_description
      << " differs from the batch fit in squared_error.";
    EXPECT_EQ(online_components.n_obs, batch_components.n_obs)
        << stream_description << " differs from the batch fit in n_obs.";
    EXPECT_NEAR(online_parameters.mu, batch_parameters.mu, parameter_tolerance)
        << stream_description << " differs from the batch fit in mu.";
    EXPECT_NEAR(
        online_parameters.sigma, batch_parameters.sigma, parameter_tolerance
    ) << stream_description
      << " differs from the batch fit in sigma.";
  }

  /**
   * @brief Folds `series` through updateComponents, starting from the
   * components of its first observation alone.
   */
  GeneralLinearLikelihoodComponents streamFromFirstObservation(
      const GeneralLinearLikelihood& likelihood,
      const std::vector<double>& series
  ) {
    GeneralLinearLikelihoodComponents components{0.0, 0.0, 0.0, 1};
    for (std::size_t i{1}; i < series.size(); i++) {
      components = likelihood.updateComponents(
          components, series.at(i), series.at(i - 1)
      );
    }
    return components;
  }

} // namespace

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
  const double initial_observation = test_vec.back();

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
  const double initial_observation = test_vec.back();

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
  const GeneralLinearLikelihood likelihood;
  const std::vector<double> series = test_support::generalLinearSeries();
  std::vector<double> extended_series = series;
  extended_series.push_back(test_support::general_linear_next_observation);

  const GeneralLinearLikelihoodComponents online_components =
      likelihood.updateComponents(
          likelihood.calculateComponents(series),
          test_support::general_linear_next_observation, series.back()
      );

  expectOnlineEqualsBatch(
      likelihood, online_components, extended_series, "1 online update"
  );
}

/**
 * @test 20 online updates from the first observation of the 21-point series
 * equal a batch fit, so the state carries across updates and the first pair
 * adds no squared error.
 */
TEST(GeneralLinearOnlineTest, StreamFromFirstObservationEqualsBatch) {
  const GeneralLinearLikelihood likelihood;
  std::vector<double> series = test_support::generalLinearSeries();
  series.push_back(test_support::general_linear_next_observation);

  expectOnlineEqualsBatch(
      likelihood, streamFromFirstObservation(likelihood, series), series,
      "A stream from the first observation"
  );
}

/**
 * @test A stream whose first observation is 0 equals a batch fit. No
 * coefficient fits the zero-lag pair, so its squared lead 25 is counted.
 */
TEST(GeneralLinearOnlineTest, StreamStartingAtZeroEqualsBatch) {
  const GeneralLinearLikelihood likelihood;
  const std::vector<double> series{0.0, 5.0, 6.0, 6.5, 7.5};

  expectOnlineEqualsBatch(
      likelihood, streamFromFirstObservation(likelihood, series), series,
      "A stream starting at 0"
  );
}

/**
 * @test Streams of the 21-point series scaled by 1e80 and by 1e-85 equal a
 * batch fit, so the squared-error update neither overflows nor underflows
 * before its division.
 */
TEST(GeneralLinearOnlineTest, StreamAtExtremeScaleEqualsBatch) {
  const double scaled_sigma_tolerance = 1e-9;
  const GeneralLinearLikelihood likelihood;
  std::vector<double> series = test_support::generalLinearSeries();
  series.push_back(test_support::general_linear_next_observation);

  for (const double scale : {1e80, 1e-85}) {
    std::vector<double> scaled_series;
    for (const double observation : series) {
      scaled_series.push_back(observation * scale);
    }
    const GeneralLinearParameters batch_parameters =
        likelihood.calculateParameters(
            likelihood.calculateComponents(scaled_series)
        );
    const GeneralLinearParameters online_parameters =
        likelihood.calculateParameters(
            streamFromFirstObservation(likelihood, scaled_series)
        );
    EXPECT_NEAR(
        online_parameters.sigma / scale, batch_parameters.sigma / scale,
        scaled_sigma_tolerance
    ) << "A stream at scale "
      << scale << " differs from the batch fit in sigma.";
  }
}
