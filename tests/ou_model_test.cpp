#include "stochastic_models/entrypoints/ou_model.h"

#include <gtest/gtest.h>

#include <cstdlib>

#include "stochastic_models/numeric_utils/helpers.h"

/**
 * @test Tests the output of the simulateOrnsteinUhlenbeck function and asserts
 * that it is of the correct size.
 */
TEST(OuModelTest, simulateOrnsteinUhlenbeckOutputTest) {
    // Declare and initialize model and test parameters.
    double mu = 0.5;
    double alpha = 0.01;
    double sigma = 0.0067;
    double start = 0.0;
    unsigned int size = 20;
    unsigned int t = 1;

    // Simulate the Ornstein-Uhlenbeck process.
    const std::vector<double> series =
        simulateOrnsteinUhlenbeck(mu, alpha, sigma, start, size, t);

    // Assert that the series is of the correct size.
    EXPECT_EQ(series.size(), size)
        << "The size of the series returned by "
           "simulateOrnsteinUhlenbeck is not equal to the expected size.";
}

/**
 * @test Tests the output of the hittingTimeDensityOrnsteinUhlenbeck function
 * and asserts that it is near the expected value.
 *
 */
TEST(OuModelTest, hittingTimeDensityOutputTest) {
    // Declare and initialize model and test parameters.
    double alpha = 0.0045;
    double mu = 0.998;
    double sigma = 0.0038;
    double first = 1.04;
    double second = 1;
    double x = 1.02;
    double tolerance = 1e-5;

    // Calculate the hitting time density.
    const double value =
        hittingTimeDensityOrnsteinUhlenbeck(x, mu, alpha, sigma, first, second);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.43046005), tolerance)
        << "The value of the hitting time density calculated by "
           "hittingTimeDensityOrnsteinUhlenbeck is not equal to the expected "
           "value.";
}
/**
 * @test Tests the output of the ornsteinUhlenbeckMaximumLikelihood function and
 * asserts that it is near the expected value.
 *
 */
TEST(OuModelTest, ornsteinUhlenbeckMaximumLikelihoodOutputTest) {
    const float tolerance = 1e-5;
    // Generate mock vector.
    const std::vector<double> test_vec{0.5, 0.25, 0.5, 0.75, 1.5, 0.5};

    // Generate likelihood estimates.
    const std::unordered_map<std::string, const double> likelihood =
        ornsteinUhlenbeckMaximumLikelihood(test_vec);

    // Expect equality for mu value.
    EXPECT_LE(abs(roundToDecimals(likelihood.at("mu"), 8) - 0.58333333),
              tolerance)
        << "ornsteinUhlenbeckMaximumLikelihood not calculating correct value "
           "mu.";
    // Expect equality for alpha value.
    EXPECT_LE(abs(roundToDecimals(likelihood.at("alpha"), 8) - 1.06784063),
              tolerance)
        << "ornsteinUhlenbeckMaximumLikelihood not calculating correct value "
           "alpha.";
    // Expect equality for sigma value.
    EXPECT_LE(abs(roundToDecimals(likelihood.at("sigma"), 8) - 0.15277777),
              tolerance)
        << "ornsteinUhlenbeckMaximumLikelihood not calculating correct value "
           "sigma.";
}
