#include <gtest/gtest.h>

#include "hitting_times/hitting_time_density.h"
#include "numeric_utils/integration.h"
#include "sde/ornstein_uhlenbeck.h"
// Demonstrate some basic assertions.
TEST(IntegrationFunctionTest, OutputTest) {
    // Run integration function.
    double alpha = 5.1;
    double mu = 0.996;
    double sigma = 1.1;
    double upper = 1.05;
    double lower = 0.8;
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(mu, alpha, sigma);
    model_func fn = &integrate_hitting_time_density;
    double value = integration(fn, model, lower, upper);
    delete model;
    EXPECT_EQ(value, 0.003993143831817661)
        << "The value of the integration is not equal to the expected value.";
}
TEST(HittingTimeDensityTest, OutputTest) {
    // Assert that the hitting_time_density function is outputting the
    // correct value.
    double alpha = 0.0045;
    double mu = 0.998;
    double sigma = 0.0038;
    double first = 1.04;
    double second = 1;
    double x = 1.02;
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(mu, alpha, sigma);
    model_func fn = &integrate_hitting_time_density;
    double value = hitting_time_density(x, fn, model, first, second);
    delete model;
    EXPECT_EQ(value, 0.4304600561745478)
        << "The value of the hitting time density is not equal to the expected "
           "value.";
}
