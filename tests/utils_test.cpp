/**
 * @file utils_test.cpp
 * @brief Unit tests for the numeric_utils module.
 *
 */
#include <gtest/gtest.h>

#include <cmath>

#include "stochastic_models/hitting_times/hitting_time_density.h"
#include "stochastic_models/numeric_utils/differentiation.h"
#include "stochastic_models/numeric_utils/helpers.h"
#include "stochastic_models/numeric_utils/integration.h"
#include "stochastic_models/numeric_utils/solvers.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"
#include "stochastic_models/trading/optimal_mean_reversion.h"

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
    EXPECT_THROW(check_function_status(status, ignore_codes),
                 std::runtime_error)
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

// Example function and type used by the brentSolver test.
struct QuadraticParams {
    double a, b, c;
};
double quadratic(double x, void *params) {
    struct QuadraticParams *p = (struct QuadraticParams *)params;

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
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
    ModelFunc fn = &integrateHittingTimeDensity;

    // Adaptive integration function.
    double value = adaptiveIntegration(fn, hitting_time_kernel, lower, upper);
    delete hitting_time_kernel;
    EXPECT_EQ(value, 0.003993143831817661)
        << "Value produced by adaptiveIntegration is not equal to the "
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
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);
    void *params = new OptimalMeanReversionParams{hitting_time_kernel, x, r};
    ModelFunc fn = &funcOptimalMeanReversionF;

    // Adaptive integration function.
    double value = semiInfiniteIntegrationUpper(fn, params, lower);

    OptimalMeanReversionParams *ptr =
        static_cast<OptimalMeanReversionParams *>(params);
    delete ptr;

    EXPECT_EQ(value, 0.30603133345784983)
        << "Value produced by semiInfiniteIntegrationUpper is not equal to "
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
    void *params = nullptr;

    // Define function to differentate.
    ModelFunc fn = [](double x, void *params) -> double { return pow(x, 2); };

    // Adaptive differentiation function.
    double value = adaptiveCentralDifferentiation(fn, params, x);
    EXPECT_LT(value - 2, 1e-5)
        << "Value produced by adaptiveCentralDifferentiation is not equal to "
           "the expected value.";
}
/**
 * @brief Test that the brentSolver function produces the
 * correct output using a toy example with a quadratic function.
 *
 */
TEST(BrentSolverFunctionTest, OutputTest) {
    // Declare and initialize model and test parameters.
    void *params = new QuadraticParams{1.0, 0.0, -5.0};

    double upper = 5;
    double lower = 0;
    float tolerance = 1e-3;

    // Initialize model and define function to solve.
    ModelFunc fn = [](double x, void *params) -> double {
        return quadratic(x, params);
    };

    // Apply brent solver.
    double value = brentSolver(fn, params, lower, upper);

    QuadraticParams *ptr = static_cast<QuadraticParams *>(params);
    delete ptr;

    EXPECT_LT(value - 2.236068, tolerance)
        << "Value produced by brentSolver is not equal to "
           "the expected value.";
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
    const float tolerance = 1e-5;

    // Create core model and optimal mean reversion instances.
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(mu, alpha, sigma);

    // Calculate upperSolverBound.
    const double value = upperSolverBound(model);
    delete model;

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.649579), tolerance)
        << "Value produced by upperSolverBound is not equal to the expected "
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
    const float tolerance = 1e-5;

    // Create core model and optimal mean reversion instances.
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(mu, alpha, sigma);

    // Calculate lowerSolverBound.
    const double value = lowerSolverBound(model);
    delete model;

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.428021), tolerance)
        << "Value produced by lowerSolverBound is not equal to the expected "
           "value.";
}
