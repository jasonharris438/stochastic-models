#include "stochastic_models/likelihood/general_linear_online.h"

#include <gtest/gtest.h>

#include <cmath>

#include "stochastic_models/likelihood/general_linear_likelihood.h"
#include "stochastic_models/numeric_utils/helpers.h"

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

    // Initialise the updater with valid values.
    const GeneralLinearLikelihood likelihood;
    GeneralLinearUpdater updater(
        likelihood.calculateLeadLagInnerProduct(test_vec),
        likelihood.calculateLagSquared(test_vec), std::pow(n_obs * sigma, 2),
        initial_observation, n_obs);

    // Calculate the new mu value.
    const double actual = updater.updateMu(new_observation);

    const double expected = -0.00139819;
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

    // Initialise the updater with valid values.
    GeneralLinearUpdater updater(0.0, 0.0, std::pow(n_obs * sigma, 2),
                                 initial_observation, n_obs);

    // Calculate the new sigma value.
    const double actual =
        updater.updateSigma(new_observation, std::exp(mu) * new_observation);

    const double expected = 45.6407;
    EXPECT_LE(abs(roundToDecimals(actual, 8) - expected), tolerance)
        << "GeneralLinearUpdater updateSigma method returning invalid value.";
}
