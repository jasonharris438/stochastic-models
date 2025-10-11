#include "stochastic_models/distributions/gaussian.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <cstdlib>
#include <gtest/gtest.h>
/**
 * @test Tests the output of the GaussianDistribution::getMean method and
 * asserts that it is equal to the mu value.
 *
 */
TEST(GaussianDistributionTest, getMeanTest) {
  // Assert that getMean returns the mean of the distribution.
  double mu = 0.996;
  double sigma = 1.1;
  GaussianDistribution* model = new GaussianDistribution(mu, sigma);
  const double value = model->getMean();
  delete model;
  EXPECT_EQ(
      value, mu
  ) << "The value returned by GaussianDistribution.getMean is not the "
       "expected value.";
}
/**
 * @test Tests the output of the GaussianDistribution::getStd method and
 * asserts that it is equal to the sigma value.
 *
 */
TEST(GaussianDistributionTest, getStdTest) {
  // Assert that getStd returns the standard deviation of the distribution.
  double mu = 0.996;
  double sigma = 1.1;

  GaussianDistribution* model = new GaussianDistribution(mu, sigma);
  const double value = model->getStd();
  delete model;
  EXPECT_EQ(
      value, sigma
  ) << "The value returned by GaussianDistribution.getStd is not the "
       "expected value.";
}
/**
 * @test Tests the output of the GaussianDistribution::cdf method and asserts
 * that it is near the expected value.
 *
 */
TEST(GaussianDistributionTest, cdfTest) {
  // Assert that Cdf returns the correct density value.
  double mu = 0.996;
  double sigma = 1.1;
  double tolerance = 1e-5;
  GaussianDistribution* model = new GaussianDistribution(mu, sigma);
  const double value = model->cdf(1.2);
  delete model;
  EXPECT_LE(abs(roundToDecimals(value, 8) - 0.57356373), tolerance)
      << "The value returned by GaussianDistribution.Cdf is not the "
         "expected value.";
}
