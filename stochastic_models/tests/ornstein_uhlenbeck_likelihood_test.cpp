#include "likelihood/ornstein_uhlenbeck_likelihood.h"

#include <gtest/gtest.h>

#include <unordered_map>
// Demonstrate some basic assertions.
TEST(OrnsteinUhlenbeckLikelihoodCalculateTest, ParameterTest) {
    // Generate mock data.
    const std::vector<double> test_vec{0.5, 0.25, 0.5, 0.75, 1.5, 0.5};
    // Generate likelihood calculator and generate estimate of mu.
    OrnsteinUhlenbeckLikelihood *likelihood = new OrnsteinUhlenbeckLikelihood();
    const std::unordered_map<std::string, const double> params{
        likelihood->calculate(test_vec)};
    delete likelihood;
    // Expect equality for mu value.
    EXPECT_EQ(params.at("mu"), 0.58333333333333337)
        << "OrnsteinUhlenbeckLikelihood not calculating correct value mu.";
    // Expect equality for alpha value.
    EXPECT_EQ(params.at("alpha"), 1.0678406300013563)
        << "OrnsteinUhlenbeckLikelihood not calculating correct value alpha.";
    // Expect equality for sigma value.
    EXPECT_EQ(params.at("sigma"), 0.15277777777777765)
        << "OrnsteinUhlenbeckLikelihood not calculating correct value sigma.";
}
