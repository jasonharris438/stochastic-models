#include "stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h"
#include "stochastic_models/numeric_utils/helpers.h"

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

// Tests that the Ornstein-Uhlenbeck likelihood class produces the correct
// parameter estimates.
TEST(OrnsteinUhlenbeckLikelihoodCalculateTest, ParameterTest) {
  const float tolerance = 1e-5;
  // Generate mock data.
  const OrnsteinUhlenbeckLikelihoodComponents components = {4.0,   3.5,  4.125,
                                                            3.375, 3.25, 6};
  // Generate likelihood calculator and generate estimates.
  OrnsteinUhlenbeckLikelihood likelihood;
  const OrnsteinUhlenbeckParameters params =
      likelihood.calculateParameters(components);

  // Expect equality for mu value.
  EXPECT_NEAR(roundToDecimals(params.mu, 8), 0.89473684, tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value mu.";
  // Expect equality for alpha value.
  EXPECT_NEAR(roundToDecimals(params.alpha, 8), 0.72054615, tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value alpha.";
  // Expect equality for sigma value.

  EXPECT_NEAR(roundToDecimals(params.sigma, 8), 0.51633524, tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value sigma.";
}

// Tests the Ornstein-Uhlenbeck likelihood class returns valid updated model
// parameter values.
TEST(OrnsteinUhlenbeckLikelihoodUpdateTest, ParameterTest) {
  const float tolerance = 1e-5;

  // Generate likelihood calculator and generate estimates.
  OrnsteinUhlenbeckLikelihood likelihood;
  const OrnsteinUhlenbeckLikelihoodComponents components = {4.0,   3.5,  4.125,
                                                            3.375, 3.25, 6};
  const OrnsteinUhlenbeckLikelihoodComponents updated_components =
      likelihood.updateComponents(components, 0.75, 1.0);

  const OrnsteinUhlenbeckParameters params =
      likelihood.calculateParameters(updated_components);

  // Expect equality for mu value.
  EXPECT_NEAR(roundToDecimals(params.mu, 8), 0.82407407, tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value mu.";
  // Expect equality for alpha value.
  EXPECT_NEAR(roundToDecimals(params.alpha, 8), 0.82667857, tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value alpha.";
  // Expect equality for sigma value.

  EXPECT_NEAR(roundToDecimals(params.sigma, 8), 0.50070985, tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value sigma.";
}

// Tests that the Ornstein-Uhlenbeck likelihood class produces the correct
// parameter estimates after an update.
TEST(OrnsteinUhlenbeckLikelihoodCalculateTest, ComponentsTest) {
  const float tolerance = 1e-5;
  // Generate mock data.
  const std::vector<double> test_vec{0.5, 0.25, 0.5, 0.75, 1.5, 1.0};
  // Generate likelihood calculator and generate estimates.
  OrnsteinUhlenbeckLikelihood likelihood;
  const OrnsteinUhlenbeckLikelihoodComponents components =
      likelihood.calculateComponents(test_vec);

  // Expect equality for the lead sum component.
  EXPECT_NEAR(roundToDecimals(components.lead_sum, 2), 4.0, tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
         "lead sum component.";
  // Expect equality for the lag sum component.
  EXPECT_NEAR(roundToDecimals(components.lag_sum, 2), 3.5, tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
         "lag sum component.";
  // Expect equality for the lead sum squared component.
  EXPECT_NEAR(
      roundToDecimals(components.lead_sum_squared, 4), 4.125, tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lead "
       "sum squared component.";
  // Expect equality for the lag sum squared component.
  EXPECT_NEAR(
      roundToDecimals(components.lag_sum_squared, 4), 3.375, tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lag "
       "sum squared component.";
  // Expect equality for the lead-lag sum product component.
  EXPECT_NEAR(
      roundToDecimals(components.lead_lag_sum_product, 2), 3.25, tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
       "lead-lag sum product component.";
  // Expect equality for the number of observations component.
  EXPECT_NEAR(roundToDecimals(components.n_obs, 2), 6.0, tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
         "number of observations component.";
}

// Tests the Ornstein-Uhlenbeck likelihood class returns valid model component
// values after an update.
TEST(OrnsteinUhlenbeckLikelihoodUpdateTest, ComponentsTest) {
  const float tolerance = 1e-5;

  // Generate likelihood calculator and generate estimates.
  OrnsteinUhlenbeckLikelihood likelihood;
  const OrnsteinUhlenbeckLikelihoodComponents components = {4.0,   3.5,  4.125,
                                                            3.375, 3.25, 6};
  const OrnsteinUhlenbeckLikelihoodComponents updated_components =
      likelihood.updateComponents(components, 0.75, 1.0);

  // Expect equality for the lead sum component.
  EXPECT_NEAR(
      roundToDecimals(updated_components.lead_sum, 2), 4.75, tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lead "
       "sum component.";
  // Expect equality for the lag sum component.
  EXPECT_NEAR(
      roundToDecimals(updated_components.lag_sum, 2), 4.5, tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lag "
       "sum component.";
  // Expect equality for the lead sum squared component.
  EXPECT_NEAR(
      roundToDecimals(updated_components.lead_sum_squared, 4), 4.6875, tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lead "
       "sum squared component.";
  // Expect equality for the lag sum squared component.
  EXPECT_NEAR(
      roundToDecimals(updated_components.lag_sum_squared, 4), 4.375, tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lag "
       "sum squared component.";
  // Expect equality for the lead-lag sum product component.
  EXPECT_NEAR(
      roundToDecimals(updated_components.lead_lag_sum_product, 2), 4.0,
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
       "lead-lag sum product component.";
  // Expect equality for the number of observations component.
  EXPECT_NEAR(roundToDecimals(updated_components.n_obs, 2), 7.0, tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
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
