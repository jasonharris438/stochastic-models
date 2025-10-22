#include "stochastic_models/entrypoints/general_sde.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <cstdlib>
#include <gtest/gtest.h>

/**
 * @test Tests the output of the updateGeneralLinearSDEModel function and
 * asserts that the values are as expected.
 */
TEST(GeneralSdeModelTest, updateGeneralLinearSDEModelOutputTest) {
  const double tolerance = 1e-5;

  // Declare and initialize model and test parameters.
  const double mu = 0.5;
  const double sigma = 0.0067;
  const double new_observation = 0.52;
  const double last_observation = 0.5;
  const double lag_squared = 0.1;
  const double lead_lag_inner_product = 0.15;
  const double squared_error = 0.002;
  const uint32_t n_obs = 100;

  // Simulate the Ornstein-Uhlenbeck process.
  const std::vector<double> series = updateGeneralLinearSDEModel(
      mu, sigma, new_observation, last_observation, lag_squared,
      lead_lag_inner_product, squared_error, n_obs
  );

  EXPECT_LE(abs(roundToDecimals(series.at(0), 8) - 0.158224), tolerance)
      << "The mu value of the parameters returned by "
         "updateGeneralLinearSDEModel is not the correct value.";
  EXPECT_LE(abs(roundToDecimals(series.at(1), 8) - 0.00788256), tolerance)
      << "The sigma value of the parameters returned by "
         "updateGeneralLinearSDEModel is not the correct value.";
}
