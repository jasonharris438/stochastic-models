#include "stochastic_models/likelihood/general_linear_likelihood.h"
#include "stochastic_models/likelihood/general_linear_online.h"

#include <gtest/gtest.h>
#include <vector>

/**
 * @file
 * @brief Unit tests for the online General Linear MLE updater.
 *
 * These tests exercise single-step updates for `mu` and `sigma` using a
 * fixed time series and assert values that were computed from a reference
 * implementation (high precision offline calculation).
 */

/**
 * @test Tests that the GeneralLinearUpdater.updateMu method returns the
 * correct value.
 *
 */
TEST(GeneralLinearOnlineTest, UpdateMuTest) {
  // Create test data.
  const std::vector<double> test_vec{1094.1, 1104.1, 1107.7, 1123.6, 1115.6,
                                     1112.7, 1118.4, 1116.9, 1127.9, 1153.2,
                                     1159.6, 1153.6, 1138.3, 1124.6, 1122.6,
                                     1134.,  1132.5, 1139.8, 1133.6, 1124.5};
  const double tolerance = 1e-5;
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

  const double expected = -0.0013319415242996231;
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
  const std::vector<double> test_vec{1094.1, 1104.1, 1107.7, 1123.6, 1115.6,
                                     1112.7, 1118.4, 1116.9, 1127.9, 1153.2,
                                     1159.6, 1153.6, 1138.3, 1124.6, 1122.6,
                                     1134.,  1132.5, 1139.8, 1133.6, 1124.5};
  const double tolerance = 5e-5;
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

  const double expected = 10.216518490188054;
  EXPECT_NEAR(actual.sigma, expected, tolerance)
      << "GeneralLinearUpdater updateSigma method returning invalid value.";
}
