#include "stochastic_models/likelihood/general_linear_likelihood.h"
#include "stochastic_models/sde/general_linear.h"
#include "support/expected_values.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <vector>
/**
 * @file
 * @brief Unit tests for the GeneralLinearModel class (mean/variance helpers).
 */

/**
 * @test Tests the return value of the getMean method.
 *
 */
TEST(GeneralLinearModelTest, GetMeanTest) {
  const GeneralLinearModel model(-0.00143647, 10.4573);
  const double expected = 0.0;
  const double actual = model.getMean();
  EXPECT_EQ(expected, actual)
      << "GeneralLinearLikelihood getMean method returning invalid value.";
}

/**
 * @test Stationary variance sigma^2 / (2|mu|) exists only for mu < 0.
 */
TEST(GeneralLinearModelTest, GetUnconditionalVarianceTest) {
  const GeneralLinearModel model(-0.00143647, 10.4573);
  EXPECT_NEAR(
      model.getUnconditionalVariance(),
      expected::general_linear_model_test::get_unconditional_variance_test, 1e-9
  ) << "getUnconditionalVariance not returning sigma^2 / (2|mu|).";
}

/**
 * @test No stationary variance exists for mu >= 0; the accessor must throw.
 */
TEST(GeneralLinearModelTest, GetUnconditionalVarianceThrowsForNonNegativeMu) {
  const GeneralLinearModel positive_mu(0.05, 0.4);
  EXPECT_THROW(positive_mu.getUnconditionalVariance(), std::domain_error)
      << "getUnconditionalVariance did not throw for mu > 0.";
  const GeneralLinearModel zero_mu(0.0, 0.4);
  EXPECT_THROW(zero_mu.getUnconditionalVariance(), std::domain_error)
      << "getUnconditionalVariance did not throw for mu = 0.";
}

/**
 * @test Conditional variance over a step t is
 * sigma^2 * (exp(2 mu t) - 1) / (2 mu).
 */
TEST(GeneralLinearModelTest, GetConditionalVarianceTest) {
  const GeneralLinearModel model(-0.00143647, 10.4573);
  EXPECT_NEAR(
      model.getConditionalVariance(1.0),
      expected::general_linear_model_test::get_conditional_variance_test_step_1,
      1e-9
  ) << "getConditionalVariance wrong at t = 1.";
  EXPECT_NEAR(
      model.getConditionalVariance(2.0),
      expected::general_linear_model_test::get_conditional_variance_test_step_2,
      1e-9
  ) << "getConditionalVariance wrong at t = 2.";
}

/**
 * @test Conditional variance at mu = 0 must take the sigma^2 * t limit, and
 * the closed form must stay continuous with that limit for a mu near 0.
 */
TEST(GeneralLinearModelTest, GetConditionalVarianceZeroMuLimitTest) {
  const GeneralLinearModel zero_mu(0.0, 0.05);
  EXPECT_NEAR(
      zero_mu.getConditionalVariance(1.0),
      expected::general_linear_model_test::
          get_conditional_variance_zero_mu_limit_test_zero_mu_step_1,
      1e-12
  ) << "GeneralLinearModel conditional variance wrong at mu = 0, t = 1.";
  EXPECT_NEAR(
      zero_mu.getConditionalVariance(2.0),
      expected::general_linear_model_test::
          get_conditional_variance_zero_mu_limit_test_zero_mu_step_2,
      1e-12
  ) << "GeneralLinearModel conditional variance wrong at mu = 0, t = 2.";

  const GeneralLinearModel near_zero_mu(1.1e-12, 0.05);
  EXPECT_NEAR(
      near_zero_mu.getConditionalVariance(1.0),
      expected::general_linear_model_test::
          get_conditional_variance_zero_mu_limit_test_near_zero_mu_step_1,
      1e-12
  ) << "GeneralLinearModel conditional variance discontinuous near mu = 0.";
}

/**
 * @test Conditional variance for a nonzero mu below 1e-12 must follow the
 * closed form. Over a step of 1000 the limit sigma^2 * t is 5e-10 relative
 * below it.
 */
TEST(GeneralLinearModelTest, GetConditionalVarianceSmallMuLongStepTest) {
  const GeneralLinearModel small_mu(5e-13, 0.05);
  EXPECT_NEAR(
      small_mu.getConditionalVariance(1000.0),
      expected::general_linear_model_test::
          get_conditional_variance_small_mu_long_step_test,
      1e-12
  ) << "GeneralLinearModel conditional variance takes the mu = 0 limit for a "
       "nonzero mu.";
}

/**
 * @test Conditional variance for a subnormal mu must equal the limit
 * sigma^2 * t, because the closed form underflows there.
 */
TEST(GeneralLinearModelTest, GetConditionalVarianceSubnormalMuTest) {
  const GeneralLinearModel subnormal_mu(
      std::numeric_limits<double>::denorm_min(), 0.05
  );
  EXPECT_NEAR(
      subnormal_mu.getConditionalVariance(0.5),
      expected::general_linear_model_test::
          get_conditional_variance_subnormal_mu_test,
      1e-15
  ) << "GeneralLinearModel conditional variance underflows for a subnormal "
       "mu.";
}
/**
 * @test size = 0 must be rejected for a contract consistent with the OU model.
 */
TEST(GeneralLinearValidationTest, simulateRejectsZeroSize) {
  const GeneralLinearModel model(0.5, 0.1);
  EXPECT_THROW(model.Simulate(0.0, 0, 1), std::invalid_argument)
      << "Simulate accepted size == 0.";
}
/**
 * @test A step that is 0, negative or not finite must be rejected. The exact
 * scheme is defined for a positive finite step only.
 */
TEST(GeneralLinearValidationTest, simulateRejectsInvalidStep) {
  const GeneralLinearModel model(0.5, 0.1);
  for (const double step :
       {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::infinity()}) {
    EXPECT_THROW(model.Simulate(0.0, 10, step), std::invalid_argument)
        << "Simulate accepted t = " << step << ".";
  }
}
/**
 * @test Per-step sample variance of simulated residuals must match the
 * closed-form conditional variance for the steps 0.5, 1 and 3. Tolerance 10%
 * (~10 standard errors of a variance estimate at n = 20000).
 */
TEST(GeneralLinearModelTest, SimulatorVarianceMatchesConditionalVariance) {
  const double mu = -0.05, sigma = 0.4;
  const GeneralLinearModel model(mu, sigma);
  for (const double step : {0.5, 1.0, 3.0}) {
    const std::vector<double> path = model.Simulate(1.0, 20000, step);
    const double growth = std::exp(mu * step);
    double sum_sq = 0.0;
    for (std::size_t i = 1; i < path.size(); ++i) {
      const double residual = path[i] - path[i - 1] * growth;
      sum_sq += residual * residual;
    }
    const double sample_variance =
        sum_sq / static_cast<double>(path.size() - 1);
    const double expected = model.getConditionalVariance(step);
    EXPECT_NEAR(sample_variance, expected, 0.1 * expected)
        << "GL simulator per-step variance off at t = " << step << ".";
  }
}

/**
 * @test Simulate-then-estimate round trip at horizons 1 and 2 with n = 50000.
 * The estimator returns mu * horizon and sigma * sqrt(horizon).
 */
TEST(GeneralLinearModelTest, SimulateEstimateRoundTrip) {
  const double mu = -0.5;
  const double sigma = 0.2;
  // Standard errors at n = 50000: mu below 0.012, sigma below 0.0015. The
  // tolerances leave a margin of 4 standard errors or more, and the seed
  // fixes the path so that a failure reproduces.
  const double mu_tolerance = 0.05;
  const double sigma_tolerance = 0.006;
  const GeneralLinearModel model(mu, sigma, std::uint64_t{40});
  const GeneralLinearLikelihood likelihood;

  for (const unsigned int horizon : {1u, 2u}) {
    const std::vector<double> path = model.Simulate(1.0, 50000, horizon);
    const GeneralLinearParameters estimate =
        likelihood.calculateParameters(likelihood.calculateComponents(path));
    const double expected_sigma =
        sigma * std::sqrt(static_cast<double>(horizon));
    EXPECT_NEAR(estimate.mu, mu * horizon, mu_tolerance)
        << "GL round trip fails to recover mu * horizon at horizon " << horizon
        << ".";
    EXPECT_NEAR(estimate.sigma, expected_sigma, sigma_tolerance)
        << "GL round trip fails to recover sigma * sqrt(horizon) at horizon "
        << horizon << ".";
  }
}

/**
 * @test 2 models built with the same seed simulate the same path and a third
 * seed changes it. Exact equality is the reproducibility claim, as in
 * GaussianSeedingTest, not a numerical result.
 */
TEST(GeneralLinearModelTest, SeededSimulateIsReproducible) {
  const GeneralLinearModel first(-0.5, 0.2, std::uint64_t{7});
  const GeneralLinearModel second(-0.5, 0.2, std::uint64_t{7});
  const GeneralLinearModel other(-0.5, 0.2, std::uint64_t{8});
  const std::vector<double> path = first.Simulate(1.0, 50, 1);
  EXPECT_EQ(path, second.Simulate(1.0, 50, 1))
      << "GL models with the same seed simulate different paths.";
  EXPECT_NE(path, other.Simulate(1.0, 50, 1))
      << "GL models with different seeds simulate the same path.";
}
