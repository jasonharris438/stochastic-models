#include "stochastic_models/likelihood/general_linear_likelihood.h"
#include "stochastic_models/likelihood/general_linear_online.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <cmath>
#include <gtest/gtest.h>
#include <iostream>

// Test that the GeneralLinearUpdater.updateMu method returns the correct value.
TEST(GeneralLinearOnlineTest, UpdateMuTest) {
  // Create test data.
  const std::vector<double> test_vec{1094.1, 1104.1, 1107.7, 1123.6, 1115.6,
                                     1112.7, 1118.4, 1116.9, 1127.9, 1153.2,
                                     1159.6, 1153.6, 1138.3, 1124.6, 1122.6,
                                     1134.,  1132.5, 1139.8, 1133.6, 1124.5};
  const size_t n_obs = test_vec.size();
  const double sigma = 10.4573;
  const double mu = -0.00143647;
  const float tolerance = 1e-5;
  const double new_observation = 1125.25;
  const double initial_observation = 1124.5;
  const double sigma_state = std::pow(n_obs * sigma, 2);

  // Calculate expected value tracking components.
  const GeneralLinearLikelihood likelihood;
  const double mu_numerator = likelihood.calculateLeadLagInnerProduct(test_vec);
  const double mu_denominator = likelihood.calculateLagSquared(test_vec);

  GeneralLinearUpdater updater(
      mu_numerator, mu_denominator, sigma_state, initial_observation, n_obs
  );

  // Calculate the new mu value.
  const double actual = updater.updateMu(new_observation);

  std::cout << "Actual: " << actual << std::endl;
  const double expected = -0.00133194;
  EXPECT_LE(abs(roundToDecimals(actual, 8) - expected), tolerance)
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
  const double sigma = 10.4573;
  const double mu = -0.00143647;
  const float tolerance = 1e-5;
  const double new_observation = 1125.25;
  const double initial_observation = 1124.5;
  const double sigma_state = std::pow(n_obs * sigma, 2);

  // Calculate expected value tracking components.
  const GeneralLinearLikelihood likelihood;
  const double mu_numerator = likelihood.calculateLeadLagInnerProduct(test_vec);
  const double mu_denominator = likelihood.calculateLagSquared(test_vec);

  // Initialise the updater with valid values.
  GeneralLinearUpdater updater(
      mu_numerator, mu_denominator, sigma_state, initial_observation, n_obs
  );

  // Calculate the new sigma value.
  const double actual = updater.updateSigma(new_observation);

  std::cout << "Actual: " << actual << std::endl;
  const double expected = 45.6394;
  EXPECT_LE(abs(roundToDecimals(actual, 8) - expected), tolerance)
      << "GeneralLinearUpdater updateSigma method returning invalid value.";
}
