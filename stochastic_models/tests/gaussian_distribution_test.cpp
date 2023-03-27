#include <gtest/gtest.h>

#include "distributions/gaussian.h"
// Demonstrate some basic assertions.
TEST(GaussianDistributionTest, getMeanTest) {
    // Assert that getMean returns the mean of the distribution.
    double mu = 0.996;
    double sigma = 1.1;
    GaussianDistribution* model = new GaussianDistribution(mu, sigma);
    const double value = model->getMean();
    delete model;
    EXPECT_EQ(value, mu)
        << "The value returned by GaussianDistribution.getMean is not the "
           "expected value.";
}
TEST(GaussianDistributionTest, getStdTest) {
    // Assert that getStd returns the standard deviation of the distribution.
    double mu = 0.996;
    double sigma = 1.1;
    GaussianDistribution* model = new GaussianDistribution(mu, sigma);
    const double value = model->getStd();
    delete model;
    EXPECT_EQ(value, sigma)
        << "The value returned by GaussianDistribution.getStd is not the "
           "expected value.";
}
TEST(GaussianDistributionTest, CdfTest) {
    // Assert that Cdf returns the correct density value.
    double mu = 0.996;
    double sigma = 1.1;
    GaussianDistribution* model = new GaussianDistribution(mu, sigma);
    const double value = model->Cdf(1.2);
    delete model;
    EXPECT_EQ(value, 0.57356373220682144)
        << "The value returned by GaussianDistribution.Cdf is not the "
           "expected value.";
}
