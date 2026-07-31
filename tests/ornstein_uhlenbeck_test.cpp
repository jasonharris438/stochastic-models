#include "stochastic_models/hitting_times/hitting_time_ornstein_uhlenbeck.h"
#include "stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h"
#include "stochastic_models/numeric_utils/helpers.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"

#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <stdexcept>
#include <vector>
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
 * @test Conditional variance over a step t must equal
 * sigma^2 * (1 - exp(-2 alpha t)) / (2 alpha).
 */
TEST(OrnsteinUhlenbeckModelTest, getConditionalVarianceOutputTest) {
  OrnsteinUhlenbeckModel model(0.5, 0.02, 0.05);
  EXPECT_NEAR(model.getConditionalVariance(1.0), 0.0024506601, 1e-9)
      << "OrnsteinUhlenbeckModel conditional variance wrong at t = 1.";
  EXPECT_NEAR(model.getConditionalVariance(2.0), 0.0048052284, 1e-9)
      << "OrnsteinUhlenbeckModel conditional variance wrong at t = 2.";
}
/**
 * @test Conditional variance at alpha = 0 must take the sigma^2 * t limit
 * branch, and must agree with that limit near the alpha = 1e-12 guard
 * threshold.
 */
TEST(OrnsteinUhlenbeckModelTest, getConditionalVarianceZeroAlphaLimitTest) {
  OrnsteinUhlenbeckModel zero_alpha(0.5, 0.0, 0.05);
  EXPECT_NEAR(zero_alpha.getConditionalVariance(1.0), 0.0025, 1e-12)
      << "OrnsteinUhlenbeckModel conditional variance wrong at alpha = 0, "
         "t = 1.";
  EXPECT_NEAR(zero_alpha.getConditionalVariance(2.0), 0.005, 1e-12)
      << "OrnsteinUhlenbeckModel conditional variance wrong at alpha = 0, "
         "t = 2.";

  // Just above the guard, so the closed-form branch is exercised at the
  // smallest alpha it ever sees.
  OrnsteinUhlenbeckModel near_zero_alpha(0.5, 1.1e-12, 0.05);
  EXPECT_NEAR(near_zero_alpha.getConditionalVariance(1.0), 0.0025, 1e-12)
      << "OrnsteinUhlenbeckModel conditional variance discontinuous across "
         "the alpha = 1e-12 guard threshold.";
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

  // sigma == 0: the exact-scheme diffusion sd is zero regardless of draws, so
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
/**
 * @test Per-step sample variance of simulated residuals must match the
 * closed-form conditional variance for t = 1 and t = 3. Tolerance 10%
 * (~10 standard errors of a variance estimate at n = 20000).
 */
TEST(OrnsteinUhlenbeckModelTest, SimulatorVarianceMatchesConditionalVariance) {
  const double mu = 1.0, alpha = 0.5, sigma = 0.3;
  OrnsteinUhlenbeckModel model(mu, alpha, sigma);
  for (const unsigned int t : {1u, 3u}) {
    const std::vector<double> path = model.Simulate(mu, 20000, t);
    const double decay = std::exp(-alpha * static_cast<double>(t));
    double sum_sq = 0.0;
    for (std::size_t i = 1; i < path.size(); ++i) {
      const double residual =
          path[i] - (path[i - 1] * decay + mu * (1 - decay));
      sum_sq += residual * residual;
    }
    const double sample_variance =
        sum_sq / static_cast<double>(path.size() - 1);
    const double expected = model.getConditionalVariance(t);
    EXPECT_NEAR(sample_variance, expected, 0.1 * expected)
        << "OU simulator per-step variance off at t = " << t << ".";
  }
}
/**
 * @test Simulate with known parameters, estimate with the MLE, recover the
 * inputs. For t = 2 the unit-step estimator sees an exact OU with
 * alpha' = 2 alpha and sigma' = sigma * sqrt(2) (mu unchanged) — see
 * docs/derivations/sde-mle-derivations.md. Tolerances >= 6 standard errors
 * at n = 50000.
 */
TEST(OrnsteinUhlenbeckModelTest, SimulateEstimateRoundTrip) {
  const double mu = 1.0, alpha = 0.7, sigma = 0.5;
  OrnsteinUhlenbeckModel model(mu, alpha, sigma);
  OrnsteinUhlenbeckLikelihood likelihood;

  const std::vector<double> path1 = model.Simulate(mu, 50000, 1);
  const OrnsteinUhlenbeckParameters p1 =
      likelihood.calculateParameters(likelihood.calculateComponents(path1));
  EXPECT_NEAR(p1.mu, mu, 0.05) << "OU round trip fails to recover mu at t=1.";
  EXPECT_NEAR(p1.alpha, alpha, 0.08)
      << "OU round trip fails to recover alpha at t=1.";
  EXPECT_NEAR(p1.sigma, sigma, 0.03)
      << "OU round trip fails to recover sigma at t=1.";

  const std::vector<double> path2 = model.Simulate(mu, 50000, 2);
  const OrnsteinUhlenbeckParameters p2 =
      likelihood.calculateParameters(likelihood.calculateComponents(path2));
  EXPECT_NEAR(p2.mu, mu, 0.05) << "OU round trip fails to recover mu at t=2.";
  EXPECT_NEAR(p2.alpha, 2 * alpha, 0.12)
      << "OU round trip fails to recover alpha * t at t=2.";
  EXPECT_NEAR(p2.sigma, sigma * std::sqrt(2.0), 0.03)
      << "OU round trip fails to recover sigma * sqrt(t) at t=2.";
}
