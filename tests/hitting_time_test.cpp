#include "stochastic_models/hitting_times/hitting_time_density.h"
#include "stochastic_models/hitting_times/hitting_time_ornstein_uhlenbeck.h"
#include "stochastic_models/numeric_utils/helpers.h"
#include "stochastic_models/numeric_utils/integration.h"

#include <cstdlib>
#include <gtest/gtest.h>
/**
 * @test Tests the output of the hittingTimeDensity function and asserts that
 * it is near the expected value.
 *
 */
TEST(HittingTimeDensityTest, hittingTimeDensityOutputTest) {
  // Declare and initialize model and test parameters.
  double alpha = 0.0045;
  double mu = 0.998;
  double sigma = 0.0038;
  double first = 1.04;
  double second = 1;
  double x = 1.02;
  double tolerance = 1e-5;
  // Create core model instance and declare function to use.
  void* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  ModelFunc fn = &integrateHittingTimeDensity;

  // Calculate the hitting time density.
  const double value =
      hittingTimeDensity(x, fn, hitting_time_kernel, first, second);
  HittingTimeOrnsteinUhlenbeck* ptr =
      static_cast<HittingTimeOrnsteinUhlenbeck*>(hitting_time_kernel);
  delete ptr;
  ptr = nullptr;
  hitting_time_kernel = nullptr;

  // Assert that the value is near the expected value.
  EXPECT_LE(abs(roundToDecimals(value, 8) - 0.43046005), tolerance)
      << "The value of the hitting time density is not equal to the expected "
         "value.";
}

/**
 * @test Tests that the hittingTimeDensity function throws the correct exception
 * when an incorrectly specified void pointer is provided.
 *
 */
TEST(HittingTimeDensityTest, CatchErrorTest) {
  double first = 1.04;
  double second = 1;
  double x = 1.02;
  void* hitting_time_kernel = new std::string{"This is not a model"};
  ModelFunc fn = &integrateHittingTimeDensity;
  EXPECT_THROW(
      {
        try {
          hittingTimeDensity(x, fn, hitting_time_kernel, first, second);
        } catch (std::runtime_error& e) {
          std::string* ptr = static_cast<std::string*>(hitting_time_kernel);
          delete ptr;
          ptr = nullptr;
          hitting_time_kernel = nullptr;
          throw;
        }
      },
      std::runtime_error
  ) << "hittingTimeDensity function must throw a runtime_error when a void "
       "pointer to the wrong type is provided.";
}
