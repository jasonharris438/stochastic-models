#include "stochastic_models/hitting_times/hitting_time_ornstein_uhlenbeck.h"
#include "stochastic_models/numeric_utils/helpers.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"

#include <gtest/gtest.h>

#include <stdexcept>
/**
 * @test Tests the output of the
 * OrnsteinUhlenbeckModel::getUnconditionalVariance method and asserts that it
 * is near the expected value.
 *
 */
TEST(OrnsteinUhlenbeckModelTest, getUnconditionalVarianceOutputTest) {
  const float tolerance = 1e-5;
  // Instantiate model and generate output of hitting time core function.
  OrnsteinUhlenbeckModel model(0.5, 0.02, 0.05);
  const double output = model.getUnconditionalVariance();

  // Assert that the value is near the expected value.
  EXPECT_NEAR(roundToDecimals(output, 8), 0.0625, tolerance)
      << "OrnsteinUhlenbeckModel not calculating "
         "correct value for getUnconditionalVariance method.";
}
/**
 * @test Tests the output of the OrnsteinUhlenbeckModel::getMean
 * method and asserts that it is near the expected value.
 *
 */
TEST(OrnsteinUhlenbeckModelTest, getMeanOutputTest) {
  const float tolerance = 1e-5;
  // Instantiate model and generate output of hitting time core function.
  OrnsteinUhlenbeckModel model(0.5, 0.02, 0.05);
  const double output = model.getMean();

  // Assert that the value is near the expected value.
  EXPECT_NEAR(roundToDecimals(output, 8), 0.5, tolerance)
      << "OrnsteinUhlenbeckModel not calculating "
         "correct value for getMean method.";
}
/**
 * @test Tests the output of the
 * HittingTimeOrnsteinUhlenbeck::hittingTimeDensityCore method and asserts that
 * it is near the expected value.
 *
 */
TEST(HittingTimeOrnsteinUhlenbeckTest, hittingTimeCoreOutputTest) {
  const float tolerance = 1e-5;
  // Instantiate model and generate output of hitting time core function.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(0.5, 0.02, 0.05);
  const double output = hitting_time_kernel.hittingTimeDensityCore(0.3);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(roundToDecimals(output, 8), 0.18637397, tolerance)
      << "OrnsteinUhlenbeckModel not calculating "
         "correct value for hitting time density.";
}
/**
 * @test Tests the output of the
 * HittingTimeOrnsteinUhlenbeck::optimalTradingFCore method and asserts that it
 * is near the expected value.
 *
 */
TEST(HittingTimeOrnsteinUhlenbeckTest, optimalTradingFCoreOutputTest) {
  const float tolerance = 1e-5;
  // Instantiate model and generate output of first passage time core
  // function F.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(0.5, 0.02, 0.05);
  const double output = hitting_time_kernel.optimalTradingFCore(0.3, 0.1, 0.02);
  // Expect equality for output.
  EXPECT_NEAR(roundToDecimals(output, 8), 0.91851228, tolerance)
      << "HittingTimeOrnsteinUhlenbeck not calculating "
         "correct value for first passage time F function.";
}
/**
 * @test Tests the output of the
 * HittingTimeOrnsteinUhlenbeck::optimalTradingGCore method and asserts that it
 * is near the expected value.
 *
 */
TEST(HittingTimeOrnsteinUhlenbeckTest, optimalTradingGCoreOutputTest) {
  const float tolerance = 1e-5;
  // Instantiate model and generate output of first passage time core
  // function G.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(0.5, 0.02, 0.05);
  const double output = hitting_time_kernel.optimalTradingGCore(0.3, 0.1, 0.02);
  // Expect equality for output.
  EXPECT_NEAR(roundToDecimals(output, 8), 1.07788415, tolerance)
      << "HittingTimeOrnsteinUhlenbeck not calculating "
         "correct value for first passage time G function.";
}
/**
 * @test Tests the output of the
 * HittingTimeOrnsteinUhlenbeck::optimalTradingLCore method and asserts that it
 * is near the expected value.
 *
 */
TEST(HittingTimeOrnsteinUhlenbeckTest, optimalTradingLCoreOutputTest) {
  const float tolerance = 1e-5;
  // Instantiate model and generate output of core function L.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(0.3, 8, 0.3);
  const double output = hitting_time_kernel.optimalTradingLCore(0.05, 0.02);
  // Expect equality for output.
  EXPECT_NEAR(roundToDecimals(output, 8), 0.298261, tolerance)
      << "HittingTimeOrnsteinUhlenbeck not calculating correct value for L "
         "function.";
}
/**
 * @test A path simulated with sigma > 0 must carry continuous Gaussian noise of
 * a plausible magnitude. A simulator must keep per-step variance near sigma^2
 * scale, so we assert the variance sits in a plausible band [1e-6, 100): the
 * lower bound rejects a dead/zero-noise path, the upper bound rejects a
 * truncation/UB explosion. Independent of the (internally random) seed.
 */
TEST(OrnsteinUhlenbeckModelTest, SimulateCarriesStochasticNoise) {
  const double start = 0.5;
  const unsigned int size = 500;
  const unsigned int t = 1;

  // sigma == 0: noise term (t * sigma * noise) is zero regardless of draws, so
  // this path is fully deterministic and identical run-to-run.
  OrnsteinUhlenbeckModel drift_only(0.5, 0.1, 0.0);
  const std::vector<double> baseline = drift_only.Simulate(start, size, t);

  // sigma > 0: a correct simulator injects continuous Gaussian noise, so the
  // path's increments must vary (positive sample variance) and the path must
  // differ from the deterministic baseline.
  OrnsteinUhlenbeckModel with_noise(0.5, 0.1, 0.3);
  const std::vector<double> noisy = with_noise.Simulate(start, size, t);

  ASSERT_EQ(baseline.size(), noisy.size());

  // Sample variance of the step-to-step increments of the noisy path.
  double mean_inc = 0.0;
  for (std::size_t i = 1; i < noisy.size(); ++i)
    mean_inc += noisy[i] - noisy[i - 1];
  mean_inc /= static_cast<double>(noisy.size() - 1);
  double var_inc = 0.0;
  for (std::size_t i = 1; i < noisy.size(); ++i) {
    const double d = (noisy[i] - noisy[i - 1]) - mean_inc;
    var_inc += d * d;
  }
  var_inc /= static_cast<double>(noisy.size() - 1);

  EXPECT_GT(
      var_inc, 1e-6
  ) << "OU Simulate increments show no stochastic variation.";
  EXPECT_LT(var_inc, 100.0)
      << "OU increments are implausibly large — integer-truncation/UB bug "
         "inflates the noise term (legitimate per-step variance ~0.1).";
  EXPECT_NE(noisy, baseline)
      << "OU noisy path is identical to the deterministic drift-only path.";
}
/**
 * @test size = 0 previously wrapped to 4294967295 in `sample(size - 1)`,
 * attempting a ~34 GB allocation. It must be rejected up front.
 */
TEST(OrnsteinUhlenbeckValidationTest, simulateRejectsZeroSize) {
  const OrnsteinUhlenbeckModel model(0.5, 0.01, 0.0067);
  EXPECT_THROW(model.Simulate(0.0, 0, 1), std::invalid_argument)
      << "Simulate accepted size == 0.";
}
