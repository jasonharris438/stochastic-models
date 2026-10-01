#include "stochastic_models/entrypoints/kca_filter.h"
#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/kalman_filter/kca.h"
#include "stochastic_models/kalman_filter/states_exceptions.h"
#include "support/assertions.h"
#include "support/expected_values.h"

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

constexpr double state_tolerance = 1e-12;

using test_support::expectMatrixNear;
using test_support::expectVectorNear;

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

  const nlohmann::json state = nlohmann::json::parse(
      getInitializedKcaState(data_series, h, q, system_dimension)
  );

  namespace initial = expected::kca_test;
  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      state.at("transition_matrix"),
      initial::get_initialized_kca_state_test_transition_matrix,
      state_tolerance, "The initialised state has the wrong transition matrix."
  ));
  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      state.at("transition_covariance"),
      initial::get_initialized_kca_state_test_transition_covariance,
      state_tolerance,
      "The initialised state has the wrong transition covariance."
  ));
  ASSERT_NO_FATAL_FAILURE(expectVectorNear(
      state.at("current_state_mean"),
      initial::get_initialized_kca_state_test_current_state_mean,
      state_tolerance, "The initialised state has the wrong current state mean."
  ));
  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      state.at("current_state_covariance"),
      initial::get_initialized_kca_state_test_current_state_covariance,
      state_tolerance,
      "The initialised state has the wrong current state covariance."
  ));
  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      state.at("observation_matrix"),
      initial::get_initialized_kca_state_test_observation_matrix,
      state_tolerance, "The initialised state has the wrong observation matrix."
  ));
  EXPECT_NEAR(
      state.at("observation_offset").get<double>(),
      initial::get_initialized_kca_state_test_observation_offset,
      state_tolerance
  ) << "The initialised state has the wrong observation offset.";
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

  const nlohmann::json updated_state = nlohmann::json::parse(
      getUpdatedKcaState(state, system_dimension, observation, innovation_sigma)
  );

  namespace updated = expected::kca_test;
  ASSERT_NO_FATAL_FAILURE(expectVectorNear(
      updated_state.at("current_state_mean"),
      updated::get_updated_kca_state_test_current_state_mean, state_tolerance,
      "The updated state has the wrong current state mean."
  ));
  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      updated_state.at("current_state_covariance"),
      updated::get_updated_kca_state_test_current_state_covariance,
      state_tolerance,
      "The updated state has the wrong current state covariance."
  ));
  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      updated_state.at("transition_matrix"),
      updated::get_updated_kca_state_test_transition_matrix, state_tolerance,
      "The updated state has the wrong transition matrix."
  ));
  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      updated_state.at("transition_covariance"),
      updated::get_updated_kca_state_test_transition_covariance,
      state_tolerance, "The updated state has the wrong transition covariance."
  ));
  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      updated_state.at("observation_matrix"),
      updated::get_updated_kca_state_test_observation_matrix, state_tolerance,
      "The updated state has the wrong observation matrix."
  ));
  EXPECT_NEAR(
      updated_state.at("observation_offset").get<double>(),
      updated::get_updated_kca_state_test_observation_offset, state_tolerance
  ) << "The updated state has the wrong observation offset.";
}

/**
 * @test A consistent dimension set that is not the fixed 3-state KCA
 * scheme must be rejected by the initialise entrypoint with the typed
 * dimensions exception.
 */
TEST(KcaValidationTest, getInitializedKcaStateRejectsNonSchemeDimensions) {
  const std::vector<double> data_series{10.5, 10.6, 10.7};
  const std::string four_state_dimensions = R"({
      "observation_covariance_columns":1,"observation_covariance_rows":1,
      "observation_matrix_columns":4,"observation_matrix_rows":1,
      "observation_offset":0.0,"state_covariance_columns":4,
      "state_covariance_rows":4,"state_mean_dimension":4})";
  EXPECT_THROW(
      getInitializedKcaState(data_series, 1.0, 0.001, four_state_dimensions),
      invalid_filter_dimensions
  ) << "getInitializedKcaState accepted a 4-state dimension set.";
}

/**
 * @test The update entrypoint must reject a non-scheme dimension set with the
 * typed dimensions exception before it reads any state.
 */
TEST(KcaValidationTest, getUpdatedKcaStateRejectsNonSchemeDimensions) {
  const std::string state = "{}";
  const std::string five_state_dimensions = R"({
      "observation_covariance_columns":1,"observation_covariance_rows":1,
      "observation_matrix_columns":5,"observation_matrix_rows":1,
      "observation_offset":0.0,"state_covariance_columns":5,
      "state_covariance_rows":5,"state_mean_dimension":5})";
  EXPECT_THROW(
      getUpdatedKcaState(state, five_state_dimensions, 10.3, 0.1),
      invalid_filter_dimensions
  ) << "getUpdatedKcaState accepted a 5-state dimension set.";
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

/**
 * @test The facade forwards an empty series to setInitialState, so it must
 * surface the same typed exception.
 */
TEST(KcaValidationTest, initialiseFilterRejectsEmptySeries) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KineticComponents kinetic_components(dimensions);
  EXPECT_THROW(
      kinetic_components.initialiseFilter({}, 1.0, 0.001),
      InvalidNumberObservationsError
  ) << "initialiseFilter accepted an empty data series.";
}

/**
 * @test The C entry point must reject an empty series with the typed
 * exception before it touches the filter state.
 */
TEST(KcaValidationTest, getInitializedKcaStateRejectsEmptySeries) {
  const std::string scheme_dimensions = R"({
      "observation_covariance_columns":1,"observation_covariance_rows":1,
      "observation_matrix_columns":3,"observation_matrix_rows":1,
      "observation_offset":0.0,"state_covariance_columns":3,
      "state_covariance_rows":3,"state_mean_dimension":3})";
  EXPECT_THROW(
      getInitializedKcaState({}, 1.0, 0.001, scheme_dimensions),
      InvalidNumberObservationsError
  ) << "getInitializedKcaState accepted an empty data series.";
}
