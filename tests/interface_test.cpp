#include <gtest/gtest.h>

#include <cstdlib>

#include "stochastic_models/interface/stochastic_models.h"
#include "stochastic_models/numeric_utils/helpers.h"
/**
 * @test Tests the output of the optimalEntryLevelLower with a stop loss
 * provided function and asserts that it is near the expected value.
 *
 */
TEST(InterfaceTest, optimalEntryLowerStopLossOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double stop_loss = 0.05;
    const double d_star = 0.136755;
    const double b_star = 0.450895;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-5;

    // Calculate a*.
    const double value = optimalEntryLevelLower(d_star, b_star, mu, alpha,
                                                sigma, stop_loss, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.118451), tolerance)
        << "Value produced by optimalEntryLevelLower function "
           "with a stop loss is not equal to the expected value.";
}
/**
 * @test Tests that the optimalEntryLevel function errors due to no solution
 * being found with the provided parameters.
 */
TEST(InterfaceTest, optimalEntryLevelNoSolutionTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 0.000116;
    const double mu = 1.818978;
    const double sigma = 0.006623;
    const double d_star = 0.436755;
    const double b_star = 0.750895;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-5;

    // Assert that the method is not implemented.
    ASSERT_THROW(optimalEntryLevel(d_star, b_star, mu, alpha, sigma, r, c),
                 std::runtime_error)
        << "optimalEntryLevel did not throw "
           "with invalid parameters.";
}
/**
 * @test Tests the output of the optimalEntryLevel with a stop loss provided
 * function and asserts that it is near the expected value.
 *
 */
TEST(InterfaceTest, optimalEntryLevelStopLossOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double stop_loss = 0.05;
    const double b_star = 0.450895;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-5;

    // Calculate d*.
    const double value =
        optimalEntryLevel(b_star, mu, alpha, sigma, stop_loss, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.136755), tolerance)
        << "Value produced by optimalEntryLevel function "
           "with a stop loss is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalExitLevel with a stop loss provided
 * function and asserts that it is near the expected value.
 *
 */
TEST(InterfaceTest, optimalExitLevelStopLossOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double stop_loss = 0.05;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-5;

    // Calculate b*.
    const double value = optimalExitLevel(mu, alpha, sigma, stop_loss, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.46912012), tolerance)
        << "Value produced by optimalExitLevel function "
           "with a stop loss is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalExitLevel function and asserts
 * that it is near the expected value.
 *
 */
TEST(InterfaceTest, optimalExitLevelOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-5;

    // Calculate b*.
    const double value = optimalExitLevel(mu, alpha, sigma, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.466836), tolerance)
        << "Value produced by optimalExitLevel function "
           "is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalExitLevelExponential function and
 * asserts that it is near the expected value.
 *
 */
TEST(InterfaceTest, optimalExitLevelExponentialOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 5;
    const double mu = 1.3499;
    const double sigma = 0.15;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Calculate b*.
    const double value = optimalExitLevelExponential(mu, alpha, sigma, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 1.4093), tolerance)
        << "Value produced by optimalExitLevelExponential function "
           "is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalEntryLevelExponential function and
 * asserts that it is near the expected value.
 *
 */
TEST(InterfaceTest, optimalEntryLevelExponentialOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 5;
    const double mu = 1.3499;
    const double sigma = 0.15;
    const double b_star = 1.4093;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Calculate d*.
    const double value =
        optimalEntryLevelExponential(b_star, mu, alpha, sigma, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 1.24096), tolerance)
        << "Value produced by optimalEntryLevelExponential function "
           "is not equal to the expected value.";
}
/**
 * @test Tests that the optimalExitLevelExponential function errors due to
 * no solution being found with the provided parameters.
 *
 */
TEST(InterfaceTest, optimalExitLevelExponentialNoSolutionTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 0.000116;
    const double mu = 1.818978;
    const double sigma = 0.006623;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Assert that the method is not implemented.
    ASSERT_THROW(optimalExitLevelExponential(mu, alpha, sigma, r, c),
                 std::runtime_error)
        << "optimalExitLevelExponential did not throw a runtime error when "
           "provided "
           "with invalid parameters.";
}
/**
 * @test Tests the output of the optimalEntryLevelLowerExponential function and
 * asserts that it is near the expected value.
 *
 */
TEST(InterfaceTest, optimalEntryLevelLowerExponentialOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 5;
    const double mu = 1.3499;
    const double sigma = 0.15;
    const double d_star = 1.24096;
    const double b_star = 1.4093;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Calculate d*.
    const double value = optimalEntryLevelLowerExponential(d_star, b_star, mu,
                                                           alpha, sigma, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 1.16016), tolerance)
        << "Value produced by optimalEntryLevelLowerExponential function "
           "is not equal to the expected value.";
}
/**
 * @test Tests the output of the optimalEntryLevel function and asserts
 * that it is near the expected value.
 *
 */
TEST(InterfaceTest, optimalEntryLevelOutputTest) {
    // Declare and initialize model and test parameters.
    const double alpha = 8;
    const double mu = 0.3;
    const double sigma = 0.3;
    const double b_star = 0.466836;
    const double c = 0.02;
    const double r = 0.05;
    const float tolerance = 1e-4;

    // Calculate d*.
    const double value = optimalEntryLevel(b_star, mu, alpha, sigma, r, c);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.116948), tolerance)
        << "Value produced by optimalEntryLevel function "
           "is not equal to the expected value.";
}
/**
 * @test Tests the output of the hittingTimeDensityOrnsteinUhlenbeck function
 * and asserts that it is near the expected value.
 *
 */
TEST(InterfaceTest, hittingTimeDensityOutputTest) {
    // Declare and initialize model and test parameters.
    double alpha = 0.0045;
    double mu = 0.998;
    double sigma = 0.0038;
    double first = 1.04;
    double second = 1;
    double x = 1.02;
    double tolerance = 1e-5;

    // Calculate the hitting time density.
    const double value =
        hittingTimeDensityOrnsteinUhlenbeck(x, mu, alpha, sigma, first, second);

    // Assert that the value is near the expected value.
    EXPECT_LE(abs(roundToDecimals(value, 8) - 0.43046005), tolerance)
        << "The value of the hitting time density calculated by "
           "hittingTimeDensityOrnsteinUhlenbeck is not equal to the expected "
           "value.";
}
/**
 * @test Tests the output of the ornsteinUhlenbeckMaximumLikelihood function and
 * asserts that it is near the expected value.
 *
 */
TEST(InterfaceTest, ornsteinUhlenbeckMaximumLikelihoodOutputTest) {
    const float tolerance = 1e-5;
    // Generate mock vector.
    const std::vector<double> test_vec{0.5, 0.25, 0.5, 0.75, 1.5, 0.5};

    // Generate likelihood estimates.
    const std::unordered_map<std::string, const double> likelihood =
        ornsteinUhlenbeckMaximumLikelihood(test_vec);

    // Expect equality for mu value.
    EXPECT_LE(abs(roundToDecimals(likelihood.at("mu"), 8) - 0.58333333),
              tolerance)
        << "ornsteinUhlenbeckMaximumLikelihood not calculating correct value "
           "mu.";
    // Expect equality for alpha value.
    EXPECT_LE(abs(roundToDecimals(likelihood.at("alpha"), 8) - 1.06784063),
              tolerance)
        << "ornsteinUhlenbeckMaximumLikelihood not calculating correct value "
           "alpha.";
    // Expect equality for sigma value.
    EXPECT_LE(abs(roundToDecimals(likelihood.at("sigma"), 8) - 0.15277777),
              tolerance)
        << "ornsteinUhlenbeckMaximumLikelihood not calculating correct value "
           "sigma.";
}
/**
 * @test Tests that the getInitializedKcaState function correctly initialises
 * a KcaStates object and returns its internal state as a JSON string with the
 * correct values.
 *
 */
TEST(InterfaceTest, getInitializedKcaStateTest) {
    // Mock data to initialise the system.
    const std::vector<double> data_series{
        10.51255, 10.51985, 10.52405, 10.4656,  10.47,    10.5403,  10.4425,
        10.3087,  10.1994,  10.1839,  10.24645, 10.1795,  10.21715, 10.14995,
        10.194,   10.22505, 10.27325, 10.25095, 10.30575, 10.27645};
    const double h{1.0};
    const double q{0.001};
    const std::string system_dimension =
        "{\"observation_covariance_columns\":1,\"observation_covariance_rows\":"
        "1,\"observation_matrix_columns\":3,\"observation_matrix_rows\":1,"
        "\"observation_offset\":0.0,\"state_covariance_columns\":3,\"state_"
        "covariance_rows\":3,\"state_mean_dimension\":3}";

    // Get the initialised internal state as a JSON string from the KCA system.
    const std::string internal_state =
        getInitializedKcaState(data_series, h, q, system_dimension);
    std::cout << internal_state << std::endl;
    EXPECT_EQ(internal_state,
              "{\"current_state_covariance\":[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,"
              "0.0,0.0]],\"current_state_mean\":[10.288741828687053,0.0,0.0],"
              "\"observation_matrix\":[[1.0,0.0,0.0]],\"observation_offset\":0."
              "0,\"transition_covariance\":[[0.12695229227341848,0.0,0.0],[0.0,"
              "0.001,0.0],[0.0,0.0,0.001]],\"transition_matrix\":[[1."
              "0011961162353782,1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0]]}")
        << "The JSON string produced by the getInitializedKcaState "
           "function is incorrect.";
}
/**
 * @test Tests that the getUpdatedKcaState function correctly performs single
 * stage update to the KCA system provided and returns a JSON object
 * representing the correct internal state.
 *
 */
TEST(InterfaceTest, getUpdatedKcaStateTest) {
    // Mock data to set the current state of the system.
    const std::string state =
        "{\"current_state_covariance\":[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,"
        "0.0,0.0]],\"current_state_mean\":[10.288741828687053,0.0,0.0],"
        "\"observation_matrix\":[[1.0,0.0,0.0]],\"observation_offset\":0."
        "0,\"transition_covariance\":[[0.12695229227341848,0.0,0.0],[0.0,"
        "0.001,0.0],[0.0,0.0,0.001]],\"transition_matrix\":[[1."
        "0011961162353782,1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0]]}";
    const std::string system_dimension =
        "{\"observation_covariance_columns\":1,\"observation_covariance_rows\":"
        "1,\"observation_matrix_columns\":3,\"observation_matrix_rows\":1,"
        "\"observation_offset\":0.0,\"state_covariance_columns\":3,\"state_"
        "covariance_rows\":3,\"state_mean_dimension\":3}";
    const double observation{10.3};
    const double innovation_sigma{0.1};

    // Get the updated internal kinetic components state.
    const std::string updated_state = getUpdatedKcaState(
        state, system_dimension, observation, innovation_sigma);

    std::cout << updated_state << std::endl;
    EXPECT_EQ(
        updated_state,
        "{\"current_state_covariance\":[[0.009269818720519449,0.0,0.0],[0.0,0."
        "001,0.0],[0.0,0.0,0.001]],\"current_state_mean\":[10.3000765492722,0."
        "0,0.0],\"observation_matrix\":[[1.0,0.0,0.0]],\"observation_offset\":"
        "0.0,\"transition_covariance\":[[0.12695229227341848,0.0,0.0],[0.0,0."
        "001,0.0],[0.0,0.0,0.001]],\"transition_matrix\":[[1.0011961162353782,"
        "1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0]]}")
        << "The JSON string produced "
           "by the getUpdatedKcaState "
           "function is incorrect.";
}
