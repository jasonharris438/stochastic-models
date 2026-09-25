/**
 * @file utils_test.cpp
 * @brief Unit tests for the numeric_utils module.
 *
 */
#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/hitting_times/hitting_time_density.h"
#include "stochastic_models/numeric_utils/differentiation.h"
#include "stochastic_models/numeric_utils/helpers.h"
#include "stochastic_models/numeric_utils/integration.h"
#include "stochastic_models/numeric_utils/solvers.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"
#include "stochastic_models/trading/optimal_mean_reversion.h"
#include "support/expected_values.h"

#include <cmath>
#include <gtest/gtest.h>
#include <limits>

/**
 * @brief Test that the check_function_status function executes without throwing
 * an exception when the status code is zero.
 *
 */
TEST(NumericUtilHelperTest, CheckFunctionStatusNoThrowTest) {
  // Declare and initialize test parameters.
  int status = 0;
  std::vector<int> ignore_codes = {1, 2, 3};

  // Check function status.
  EXPECT_NO_THROW(check_function_status(status, ignore_codes))
      << "check_function_status failed to execute without throwing an "
         "exception when status == 0.";
}

/**
 * @brief Test that the check_function_status function throws an exception when
 * the status code is not zero and not ignored.
 *
 */
TEST(NumericUtilHelperTest, CheckFunctionStatusThrowTest) {
  // Declare and initialize test parameters.
  int status = 4;
  std::vector<int> ignore_codes = {1, 2, 3};

  // Check function status.
  EXPECT_THROW(check_function_status(status, ignore_codes), std::runtime_error)
      << "check_function_status failed to throw an exception when expected "
         "to.";
}

/**
 * @brief Test that the check_function_status function executes without throwing
 * an exception when the status code is ignored.
 *
 */
TEST(NumericUtilHelperTest, CheckFunctionStatusIgnoreErrorTest) {
  // Declare and initialize test parameters.
  int status = 1;
  std::vector<int> ignore_codes = {1, 2, 3};

  // Check function status.
  EXPECT_NO_THROW(check_function_status(status, ignore_codes))
      << "check_function_status threw an exception when expected "
         "to ignore error code.";
}

/**
 * @brief Test that the check_function_status function executes without throwing
 * an exception when the status code is zero and no ignore codes are provided.
 *
 */
TEST(NumericUtilHelperTest, CheckFunctionStatusNoIgnoreTest) {
  // Declare and initialize test parameters.
  int status = 0;
  std::vector<int> ignore_codes = {};

  // Check function status.
  EXPECT_NO_THROW(check_function_status(status, ignore_codes))
      << "check_function_status threw an exception when expected "
         "to complete with no ignores.";
}

/**
 * @test Tests that the check_finite_result function executes without throwing
 * when the value is finite.
 *
 */
TEST(NumericUtilHelperTest, CheckFiniteResultNoThrowTest) {
  EXPECT_NO_THROW(check_finite_result(1.5, "test_routine"))
      << "check_finite_result threw an exception for a finite value.";
}

/**
 * @test Tests that the check_finite_result function throws
 * NonFiniteResultError when the value is NaN or Inf.
 *
 */
TEST(NumericUtilHelperTest, CheckFiniteResultThrowTest) {
  EXPECT_THROW(
      check_finite_result(std::nan(""), "test_routine"), NonFiniteResultError
  ) << "check_finite_result did not throw for a NaN value.";
  EXPECT_THROW(
      check_finite_result(
          std::numeric_limits<double>::infinity(), "test_routine"
      ),
      NonFiniteResultError
  ) << "check_finite_result did not throw for an Inf value.";
}

/**
 * @test Tests that the exception message thrown by check_finite_result
 * names the routine passed to it.
 *
 */
TEST(NumericUtilHelperTest, CheckFiniteResultMessageNamesRoutineTest) {
  try {
    check_finite_result(std::nan(""), "test_routine");
    FAIL() << "check_finite_result did not throw for a NaN value.";
  } catch (const NonFiniteResultError& error) {
    EXPECT_NE(std::string(error.what()).find("test_routine"), std::string::npos)
        << "check_finite_result message does not name the routine.";
  }
}

/**
 * @test Tests that check_minimum_observations executes without throwing when
 * the series meets the minimum.
 */
TEST(NumericUtilHelperTest, CheckMinimumObservationsNoThrowTest) {
  EXPECT_NO_THROW(check_minimum_observations({1.0, 2.0}, 2, "test_routine"))
      << "check_minimum_observations threw for a series at the minimum.";
}

/**
 * @test Tests that check_minimum_observations throws
 * InvalidNumberObservationsError when the series is shorter than the minimum.
 */
TEST(NumericUtilHelperTest, CheckMinimumObservationsThrowTest) {
  EXPECT_THROW(
      check_minimum_observations({}, 1, "test_routine"),
      InvalidNumberObservationsError
  ) << "check_minimum_observations did not throw for an empty series.";
  EXPECT_THROW(
      check_minimum_observations({1.0}, 2, "test_routine"),
      InvalidNumberObservationsError
  ) << "check_minimum_observations did not throw for a series below the "
       "minimum.";
}

/**
 * @test Tests that the exception message thrown by check_minimum_observations
 * names the routine passed to it.
 */
TEST(NumericUtilHelperTest, CheckMinimumObservationsMessageNamesRoutineTest) {
  try {
    check_minimum_observations({}, 1, "test_routine");
    FAIL() << "check_minimum_observations did not throw for an empty series.";
  } catch (const InvalidNumberObservationsError& error) {
    EXPECT_NE(std::string(error.what()).find("test_routine"), std::string::npos)
        << "check_minimum_observations message does not name the routine.";
  }
}

// Example function and type used by the brentSolver test.
struct QuadraticParams {
  double a, b, c;
};
double quadratic(double x, void* params) {
  struct QuadraticParams* p = (struct QuadraticParams*)params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  return (a * x + b) * x + c;
}

/**
 * @brief Test that the adaptiveIntegration function produces the correct
 * output.
 *
 */
TEST(AdaptiveIntegrationFunctionTest, OutputTest) {
  // Declare and initialize model and test parameters.
  double alpha = 5.1;
  double mu = 0.996;
  double sigma = 1.1;
  double upper = 1.05;
  double lower = 0.8;

  // Initialize model and define function to integrate.
  HittingTimeOrnsteinUhlenbeck hitting_time_kernel(mu, alpha, sigma);
  ModelFunc fn = &integrateHittingTimeDensity;

  // Adaptive integration function.
  double value = adaptiveIntegration(fn, &hitting_time_kernel, lower, upper);
  EXPECT_NEAR(
      value, expected::adaptive_integration_function_test::output_test, 5e-9
  ) << "Value produced by adaptiveIntegration is not equal to the "
       "expected value.";
}
/**
 * @brief Test that the semiInfiniteIntegrationUpper function produces the
 * correct output.
 *
 */
TEST(SemiInfiniteIntegrationFunctionTest, OutputTest) {
  // Declare and initialize model and test parameters.
  double alpha = 5.1;
  double mu = 0.996;
  double sigma = 1.1;
  double lower = 0.8;
  double x = 0.9;
  double r = 0.03;

  // Initialize model and define function to integrate.
  auto* hitting_time_kernel =
      new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
  OptimalMeanReversionParams params =
      OptimalMeanReversionParams{hitting_time_kernel, x, r};
  ModelFunc fn = &funcOptimalMeanReversionF;

  // Adaptive integration function.
  double value = semiInfiniteIntegrationUpper(fn, &params, lower);

  EXPECT_NEAR(
      value, expected::semi_infinite_integration_function_test::output_test,
      5e-7
  ) << "Value produced by semiInfiniteIntegrationUpper is not equal to "
       "the expected value.";
}
/**
 * @brief Test that the adaptiveCentralDifferentiation function produces the
 * correct output.
 *
 */
TEST(AdaptiveCentralDifferentiationFunctionTest, OutputTest) {
  // Declare and initialize test parameters.
  double x = 1;
  void* params = nullptr;

  // Define function to differentate.
  ModelFunc fn = [](double x, void*) -> double { return pow(x, 2); };

  // Adaptive differentiation function.
  double value = adaptiveCentralDifferentiation(fn, params, x);
  EXPECT_NEAR(
      value,
      expected::adaptive_central_differentiation_function_test::output_test,
      1e-5
  ) << "Value produced by adaptiveCentralDifferentiation is not equal to "
       "the expected value.";
}
double alwaysNan(double, void*) {
  return std::nan("");
}

/**
 * @test Tests that adaptiveCentralDifferentiation throws
 * NonFiniteResultError when the function returns NaN.
 *
 */
TEST(DifferentiationValidationTest, DerivativeOfNanThrowTest) {
  double x = 1;

  ModelFunc fn = &alwaysNan;

  EXPECT_THROW(
      adaptiveCentralDifferentiation(fn, nullptr, x), NonFiniteResultError
  ) << "adaptiveCentralDifferentiation did not throw NonFiniteResultError "
       "for a NaN function.";
}

// NaN inside (1, 2); changes sign over [0, 3].
double nanInsideBracket(double x, void*) {
  if (x > 1.0 && x < 2.0) {
    return std::nan("");
  }
  return x - 1.5;
}

/**
 * @brief Test that the brentSolver function produces the
 * correct output using a toy example with a quadratic function.
 *
 */
TEST(BrentSolverFunctionTest, OutputTest) {
  // Declare and initialize model and test parameters.
  QuadraticParams params = QuadraticParams{1.0, 0.0, -5.0};

  double upper = 5;
  double lower = 0;
  const double tolerance = 3e-4;

  // Initialize model and define function to solve.
  ModelFunc fn = [](double x, void* params) -> double {
    return quadratic(x, params);
  };

  // Apply brent solver.
  double value = brentSolver(fn, &params, lower, upper);

  EXPECT_NEAR(
      value, expected::brent_solver_function_test::output_test, tolerance
  ) << "Value produced by brentSolver is not equal to "
       "the expected value.";
}
/**
 * @test Tests that brentSolver throws RootNotBracketedError when the function
 * does not change sign over the bracket.
 *
 */
TEST(BrentSolverValidationTest, BracketWithoutRootThrowTest) {
  QuadraticParams params = QuadraticParams{1.0, 0.0, 1.0};
  double lower = 0;
  double upper = 3;

  ModelFunc fn = [](double x, void* params) -> double {
    return quadratic(x, params);
  };

  EXPECT_THROW(brentSolver(fn, &params, lower, upper), RootNotBracketedError)
      << "brentSolver did not throw RootNotBracketedError for a bracket "
         "with no sign change.";
}
/**
 * @test Tests that brentSolver throws NoSolutionError when the objective
 * returns NaN at an iteration point.
 *
 */
TEST(BrentSolverValidationTest, NanObjectiveThrowTest) {
  double lower = 0;
  double upper = 3;

  ModelFunc fn = &nanInsideBracket;

  EXPECT_THROW(brentSolver(fn, nullptr, lower, upper), NoSolutionError)
      << "brentSolver did not throw NoSolutionError for a NaN objective "
         "inside the bracket.";
}

// Sign change at zero; bisection cannot meet the tolerance in 100
// iterations on a 1e23-wide bracket.
double stepFunction(double x, void*) {
  return (x < 0.0) ? -1.0 : 1.0;
}

double identityFunction(double x, void*) {
  return x;
}

/**
 * @test Tests that brentSolver throws SolverConvergenceError when the
 * interval cannot reach the tolerance within the iteration limit.
 *
 */
TEST(BrentSolverValidationTest, NonConvergenceThrowTest) {
  double lower = -1e23;
  double upper = 1e23;

  ModelFunc fn = &stepFunction;

  EXPECT_THROW(brentSolver(fn, nullptr, lower, upper), SolverConvergenceError)
      << "brentSolver did not throw SolverConvergenceError for a bracket "
         "too wide to converge.";
}

/**
 * @test Tests that brentSolver converges without throwing when the root is
 * at zero.
 *
 */
TEST(BrentSolverValidationTest, RootAtZeroConvergesTest) {
  double lower = -1;
  double upper = 2;

  ModelFunc fn = &identityFunction;

  double value = 1;
  EXPECT_NO_THROW(value = brentSolver(fn, nullptr, lower, upper))
      << "brentSolver threw an exception for a root at zero.";
  EXPECT_NEAR(value, 0.0, 1e-8)
      << "brentSolver root for f(x) = x is not near zero.";
}

// NaN inside (0.4, 0.6); the Kronrod midpoint node at 0.5 always hits it.
double nanMidInterval(double x, void*) {
  if (x > 0.4 && x < 0.6) {
    return std::nan("");
  }
  return 1.0;
}

double cosineFunction(double x, void*) {
  return std::cos(x);
}

double sineFunction(double x, void*) {
  return std::sin(x);
}

/**
 * @test Tests that adaptiveIntegration throws NonFiniteResultError when the
 * integrand returns NaN inside the interval.
 *
 */
TEST(IntegrationValidationTest, NanIntegrandThrowTest) {
  double lower = 0;
  double upper = 1;

  ModelFunc fn = &nanMidInterval;

  EXPECT_THROW(
      adaptiveIntegration(fn, nullptr, lower, upper), NonFiniteResultError
  ) << "adaptiveIntegration did not throw NonFiniteResultError for a NaN "
       "integrand.";
}

/**
 * @test Tests that adaptiveIntegration accepts a cancelling integrand whose
 * exact result is zero, where a purely relative bound would collapse.
 *
 */
TEST(IntegrationValidationTest, CancellingIntegrandAcceptedTest) {
  double lower = -1;
  double upper = 1;

  ModelFunc fn = &sineFunction;

  double value = 1;
  EXPECT_NO_THROW(value = adaptiveIntegration(fn, nullptr, lower, upper))
      << "adaptiveIntegration threw an exception for a cancelling integrand.";
  EXPECT_NEAR(value, 0.0, 1e-9)
      << "adaptiveIntegration result for a cancelling integrand is not near "
         "zero.";
}

/**
 * @test Tests that semiInfiniteIntegrationUpper throws
 * IntegrationToleranceError when the error estimate is over the accepted
 * bound.
 *
 */
TEST(IntegrationValidationTest, ToleranceThrowTest) {
  double lower = 0;

  ModelFunc fn = &cosineFunction;

  EXPECT_THROW(
      semiInfiniteIntegrationUpper(fn, nullptr, lower),
      IntegrationToleranceError
  ) << "semiInfiniteIntegrationUpper did not throw IntegrationToleranceError "
       "for a divergent integral.";
}

/**
 * @test Tests the output of the upperSolverBound function is near the expected
 * value.
 *
 */
TEST(SolverBoundsTest, upperSolverBoundOutputTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 16.6677;
  const double mu = 0.5388;
  const double sigma = 0.1599;
  const double tolerance = 1e-9;

  // Create core model and optimal mean reversion instances.
  OrnsteinUhlenbeckModel model(mu, alpha, sigma);

  // Calculate upperSolverBound.
  const double value = upperSolverBound(&model);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value, expected::solver_bounds_test::upper_solver_bound_output_test,
      tolerance
  ) << "Value produced by upperSolverBound is not equal to the expected "
       "value.";
}
/**
 * @test Tests the output of the lowerSolverBound function is near the expected
 * value.
 *
 */
TEST(SolverBoundsTest, lowerSolverBoundOutputTest) {
  // Declare and initialize model and test parameters.
  const double alpha = 16.6677;
  const double mu = 0.5388;
  const double sigma = 0.1599;
  const double tolerance = 1e-9;

  // Create core model and optimal mean reversion instances.
  OrnsteinUhlenbeckModel model(mu, alpha, sigma);

  // Calculate lowerSolverBound.
  const double value = lowerSolverBound(&model);

  // Assert that the value is near the expected value.
  EXPECT_NEAR(
      value, expected::solver_bounds_test::lower_solver_bound_output_test,
      tolerance
  ) << "Value produced by lowerSolverBound is not equal to the expected "
       "value.";
}
