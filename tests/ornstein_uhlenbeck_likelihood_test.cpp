#include "stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <cstdlib>
#include <gtest/gtest.h>
#include <unordered_map>
// Demonstrate some basic assertions.
TEST(OrnsteinUhlenbeckLikelihoodCalculateTest, ParameterTest) {
  const float tolerance = 1e-5;
  // Generate mock data.
  const std::vector<double> test_vec{0.5, 0.25, 0.5, 0.75, 1.5, 0.5};
  // Generate likelihood calculator and generate estimates.
  OrnsteinUhlenbeckLikelihood* likelihood = new OrnsteinUhlenbeckLikelihood();
  const std::unordered_map<std::string, const double> params{
      likelihood->calculate(test_vec)
  };
  delete likelihood;
  // Expect equality for mu value.
  EXPECT_LE(abs(roundToDecimals(params.at("mu"), 8) - 0.58333333), tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value mu.";
  // Expect equality for alpha value.
  EXPECT_LE(abs(roundToDecimals(params.at("alpha"), 8) - 1.06784063), tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value alpha.";
  // Expect equality for sigma value.
  EXPECT_LE(abs(roundToDecimals(params.at("sigma"), 8) - 0.15277777), tolerance)
      << "OrnsteinUhlenbeckLikelihood not calculating correct value sigma.";
}
