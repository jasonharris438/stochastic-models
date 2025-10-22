#include "stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <cstdlib>
#include <gtest/gtest.h>
#include <unordered_map>

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
  EXPECT_LE(abs(roundToDecimals(params.mu, 8) - 0.85), tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value mu.";
  // Expect equality for alpha value.
  EXPECT_LE(abs(roundToDecimals(params.alpha, 8) - 0.374693), tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value alpha.";
  // Expect equality for sigma value.

  EXPECT_LE(abs(roundToDecimals(params.sigma, 8) - 0.359877), tolerance)
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
  EXPECT_LE(abs(roundToDecimals(params.mu, 8) - 0.741667), tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value mu.";
  // Expect equality for alpha value.
  EXPECT_LE(abs(roundToDecimals(params.alpha, 8) - 0.448549), tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value alpha.";
  // Expect equality for sigma value.

  EXPECT_LE(abs(roundToDecimals(params.sigma, 8) - 0.264907), tolerance)
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
  EXPECT_LE(abs(roundToDecimals(components.lead_sum, 2) - 4.0), tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
         "lead sum component.";
  // Expect equality for the lag sum component.
  EXPECT_LE(abs(roundToDecimals(components.lag_sum, 2) - 3.5), tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
         "lag sum component.";
  // Expect equality for the lead sum squared component.
  EXPECT_LE(
      abs(roundToDecimals(components.lead_sum_squared, 4) - 4.125), tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lead "
       "sum squared component.";
  // Expect equality for the lag sum squared component.
  EXPECT_LE(
      abs(roundToDecimals(components.lag_sum_squared, 4) - 3.375), tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lag "
       "sum squared component.";
  // Expect equality for the lead-lag sum product component.
  EXPECT_LE(
      abs(roundToDecimals(components.lead_lag_sum_product, 2) - 3.25), tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
       "lead-lag sum product component.";
  // Expect equality for the number of observations component.
  EXPECT_LE(abs(roundToDecimals(components.n_obs, 2) - 6.0), tolerance)
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
  EXPECT_LE(
      abs(roundToDecimals(updated_components.lead_sum, 2) - 4.75), tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lead "
       "sum component.";
  // Expect equality for the lag sum component.
  EXPECT_LE(
      abs(roundToDecimals(updated_components.lag_sum, 2) - 4.5), tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lag "
       "sum component.";
  // Expect equality for the lead sum squared component.
  EXPECT_LE(
      abs(roundToDecimals(updated_components.lead_sum_squared, 4) - 4.6875),
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lead "
       "sum squared component.";
  // Expect equality for the lag sum squared component.
  EXPECT_LE(
      abs(roundToDecimals(updated_components.lag_sum_squared, 4) - 4.375),
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the lag "
       "sum squared component.";
  // Expect equality for the lead-lag sum product component.
  EXPECT_LE(
      abs(roundToDecimals(updated_components.lead_lag_sum_product, 2) - 4.0),
      tolerance
  ) << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
       "lead-lag sum product component.";
  // Expect equality for the number of observations component.
  EXPECT_LE(abs(roundToDecimals(updated_components.n_obs, 2) - 7.0), tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value for the "
         "number of observations component.";
}
