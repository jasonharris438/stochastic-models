#include "stochastic_models/likelihood/general_linear_likelihood.h"
#include "stochastic_models/likelihood/general_linear_online.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <cmath>
#include <gtest/gtest.h>

/**
 * @file
 * @brief Unit tests for the online General Linear MLE updater.
 *
 * These tests exercise single-step updates for `mu` and `sigma` using a
 * fixed time series and assert values that were computed from a reference
 * implementation (high precision offline calculation).
 */

// Test that the GeneralLinearUpdater.updateMu method returns the correct value.
TEST(GeneralLinearOnlineTest, UpdateMuTest) {
  // Create test data.
  const std::vector<double> test_vec{1094.1, 1104.1, 1107.7, 1123.6, 1115.6,
                                     1112.7, 1118.4, 1116.9, 1127.9, 1153.2,
                                     1159.6, 1153.6, 1138.3, 1124.6, 1122.6,
                                     1134.,  1132.5, 1139.8, 1133.6, 1124.5};
  const size_t n_obs = test_vec.size();
  const float tolerance = 1e-5;
  const double new_observation = 1125.25;
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

  const double expected = -0.00133194;
  EXPECT_LE(abs(roundToDecimals(actual.mu, 8) - expected), tolerance)
      << "GeneralLinearUpdater updateMu method returning invalid value.";
}

// Test that the GeneralLinearUpdater.updateSigma method returns the correct
// value.
TEST(GeneralLinearOnlineTest, UpdateSigmaTest) {
  // Create test data.
  const std::vector<double> test_vec{1094.1, 1104.1, 1107.7, 1123.6, 1115.6,
                                     1112.7, 1118.4, 1116.9, 1127.9, 1153.2,
                                     1159.6, 1153.6, 1138.3, 1124.6, 1122.6,
                                     1134.,  1132.5, 1139.8, 1133.6, 1124.5};
  const size_t n_obs = test_vec.size();
  const float tolerance = 1e-5;
  const double new_observation = 1125.25;
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

  const double expected = 10.2165;
  EXPECT_LE(abs(roundToDecimals(actual.sigma, 4) - expected), tolerance)
      << "GeneralLinearUpdater updateSigma method returning invalid value.";
}
