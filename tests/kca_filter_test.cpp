#include "stochastic_models/entrypoints/kca_filter.h"
#include "stochastic_models/kalman_filter/states_exceptions.h"

#include <cstdlib>
#include <gtest/gtest.h>

/**
 * @test Tests that the getInitializedKcaState function correctly initialises
 * a KcaStates object and returns its internal state as a JSON string with the
 * correct values.
 *
 */
TEST(KcaTest, getInitializedKcaStateTest) {
  // Mock data to initialise the system.
  const std::vector<double> data_series{10.51255, 10.51985, 10.52405, 10.4656,
                                        10.47,    10.5403,  10.4425,  10.3087,
                                        10.1994,  10.1839,  10.24645, 10.1795,
                                        10.21715, 10.14995, 10.194,   10.22505,
                                        10.27325, 10.25095, 10.30575, 10.27645};
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
  EXPECT_EQ(
      internal_state,
      "{\"current_state_covariance\":[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,"
      "0.0,0.0]],\"current_state_mean\":[10.27645,0.0,0.0],"
      "\"observation_matrix\":[[1.0,0.0,0.0]],\"observation_offset\":0."
      "0,\"transition_covariance\":[[5e-05,0.000125,0.00016666666666666666],"
      "[0.000125,0.0003333333333333333,0.0005],[0.00016666666666666666,0.0005,"
      "0.001]],"
      "\"transition_matrix\":[[1.0,1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0]]}"
  ) << "The JSON string produced by the getInitializedKcaState "
       "function is incorrect.";
}
/**
 * @test Tests that the getUpdatedKcaState function correctly performs single
 * stage update to the KCA system provided and returns a JSON object
 * representing the correct internal state.
 *
 */
TEST(KcaTest, getUpdatedKcaStateTest) {
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
      state, system_dimension, observation, innovation_sigma
  );

  std::cout << updated_state << std::endl;
  EXPECT_EQ(
      updated_state,
      "{\"current_state_covariance\":[[0.009269818720519449,0.0,0.0],[0.0,0."
      "001,0.0],[0.0,0.0,0.001]],\"current_state_mean\":[10.3000765492722,0."
      "0,0.0],\"observation_matrix\":[[1.0,0.0,0.0]],\"observation_offset\":"
      "0.0,\"transition_covariance\":[[0.12695229227341848,0.0,0.0],[0.0,0."
      "001,0.0],[0.0,0.0,0.001]],\"transition_matrix\":[[1.0011961162353782,"
      "1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0]]}"
  ) << "The JSON string produced "
       "by the getUpdatedKcaState "
       "function is incorrect.";
}

/**
 * @test Dimensions alone were previously sufficient for an out-of-bounds
 * write: setInitialState writes a fixed 3-state kinematic system into
 * matrices sized from caller dimensions. Non-scheme dimensions must be
 * rejected with json_parse_error.
 */
TEST(KcaValidationTest, getInitializedKcaStateRejectsNonSchemeDimensions) {
  const std::vector<double> data_series{10.5, 10.6, 10.7};
  const std::string undersized_dimensions = R"({
      "observation_covariance_columns":1,"observation_covariance_rows":1,
      "observation_matrix_columns":3,"observation_matrix_rows":1,
      "observation_offset":0.0,"state_covariance_columns":3,
      "state_covariance_rows":1,"state_mean_dimension":3})";
  EXPECT_THROW(
      getInitializedKcaState(data_series, 1.0, 0.001, undersized_dimensions),
      json_parse_error
  ) << "getInitializedKcaState accepted state_covariance_rows == 1.";
}

/**
 * @test The update entrypoint must also reject non-scheme dimensions before
 * deserializing any state.
 */
TEST(KcaValidationTest, getUpdatedKcaStateRejectsNonSchemeDimensions) {
  const std::string state = "{}";
  const std::string oversized_dimensions = R"({
      "observation_covariance_columns":1,"observation_covariance_rows":1,
      "observation_matrix_columns":3,"observation_matrix_rows":1,
      "observation_offset":0.0,"state_covariance_columns":5,
      "state_covariance_rows":5,"state_mean_dimension":5})";
  EXPECT_THROW(
      getUpdatedKcaState(state, oversized_dimensions, 10.3, 0.1),
      json_parse_error
  ) << "getUpdatedKcaState accepted 5x5 state covariance dimensions.";
}

/**
 * @test A state JSON whose arrays disagree with the (valid) dimensions must
 * surface json_parse_error through the update entrypoint.
 */
TEST(KcaValidationTest, getUpdatedKcaStateRejectsShapeMismatchedState) {
  const std::string oversized_state = R"({"current_state_covariance":
      [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]],
      "current_state_mean":[10.28,0.0,0.0],
      "observation_matrix":[[1.0,0.0,0.0]],"observation_offset":0.0,
      "transition_covariance":[[0.1,0.0,0.0],[0.0,0.001,0.0],[0.0,0.0,0.001]],
      "transition_matrix":[[1.0,1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0],
      [9.0,9.0,9.0]]})";
  const std::string scheme_dimensions = R"({
      "observation_covariance_columns":1,"observation_covariance_rows":1,
      "observation_matrix_columns":3,"observation_matrix_rows":1,
      "observation_offset":0.0,"state_covariance_columns":3,
      "state_covariance_rows":3,"state_mean_dimension":3})";
  EXPECT_THROW(
      getUpdatedKcaState(oversized_state, scheme_dimensions, 10.3, 0.1),
      json_parse_error
  ) << "getUpdatedKcaState accepted a 4-row transition_matrix.";
}
