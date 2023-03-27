#include <gtest/gtest.h>

#include "numeric_utils/integration.h"
// Demonstrate some basic assertions.
TEST(IntegrationFunctionTest, OutputTest) {
    // Run integration function.
    double alpha = 5.1;
    double mu = 0.996;
    double sigma = 1.1;
    double upper = 1.05;
    double lower = 0.8;
    double value = integration(alpha, mu, sigma, lower, upper);
    EXPECT_EQ(value, 0.003993143831817661)
        << "The value of the integration is not equal to the expected value.";
}
TEST(HittingTimeDensityTest, OutputTest) {
    // Assert that the hitting_time_density function is outputting the correct
    // value.
    double alpha = 5.1;
    double mu = 0.996;
    double sigma = 1.1;
    double first = 1.05;
    double second = 0.8;
    double x = 1.02;
    double value = hitting_time_density(x, alpha, mu, sigma, first, second);
    EXPECT_EQ(value, 0.88442954018142272)
        << "The value of the hitting time density is not equal to the expected "
           "value.";
}
