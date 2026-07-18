#include "stochastic_models/kalman_filter/adapters.h"
#include "stochastic_models/kalman_filter/states.h"
#include "stochastic_models/kalman_filter/states_exceptions.h"
#include "stochastic_models/kalman_filter/type_conversion.h"

#include <gtest/gtest.h>

/**
 * @brief Test that the FilterSystemDimensionsJsonAdapter deserialize method
 * returns the correct result.
 */
TEST(AdaptersTest, FilterSystemDimensionsJsonAdapterDeserializeTest) {
  const std::string state =
      "{\"state_mean_dimension\": 8, \"state_covariance_rows\": 7, "
      "\"state_covariance_columns\": 6, \"observation_matrix_rows\": 5, "
      "\"observation_matrix_columns\": 4, \"observation_covariance_rows\": "
      "3, \"observation_covariance_columns\": 2, \"observation_offset\": "
      "1.0}";
  FilterSystemDimensionsJsonAdapter adapter;

  const FilterSystemDimensions dimensions = adapter.deserialize(state);
  EXPECT_EQ(dimensions.state_mean_dimension, 8)
      << "The state mean dimension produced by the "
         "FilterSystemDimensions.deserialize method is incorrect.";
  EXPECT_EQ(dimensions.state_covariance_rows, 7)
      << "The state covariance rows produced by the "
         "FilterSystemDimensions.deserialize method is incorrect.";
  EXPECT_EQ(dimensions.state_covariance_columns, 6)
      << "The state covariance columns produced by the "
         "FilterSystemDimensions.deserialize method is incorrect.";
  EXPECT_EQ(dimensions.observation_matrix_rows, 5)
      << "The observation matrix rows produced by the "
         "FilterSystemDimensions.deserialize method is incorrect.";
  EXPECT_EQ(dimensions.observation_matrix_columns, 4)
      << "The observation matrix columns produced by the "
         "FilterSystemDimensions.deserialize method is incorrect.";
  EXPECT_EQ(dimensions.observation_covariance_rows, 3)
      << "The observation covariance rows produced by the "
         "FilterSystemDimensions.deserialize method is incorrect.";
  EXPECT_EQ(dimensions.observation_covariance_columns, 2)
      << "The observation covariance columns produced by the "
         "FilterSystemDimensions.deserialize method is incorrect.";
  EXPECT_EQ(dimensions.observation_offset, 1.0)
      << "The observation offset produced by the "
         "FilterSystemDimensions.deserialize method is incorrect.";
}

/**
 * @brief Test that the FilterSystemDimensionsJsonAdapter serialize method
 * returns the correct result.
 */
TEST(AdaptersTest, FilterSystemDimensionsJsonAdapterSerializeTest) {
  const FilterSystemDimensions dimensions = {8, 7, 6, 5, 4, 3, 2, 1.0};
  const std::string state =
      "{\"observation_covariance_columns\":2,\"observation_covariance_rows\":"
      "3,\"observation_matrix_columns\":4,\"observation_matrix_rows\":5,"
      "\"observation_offset\":1.0,\"state_covariance_columns\":6,\"state_"
      "covariance_rows\":7,\"state_mean_dimension\":8}";
  FilterSystemDimensionsJsonAdapter adapter;
  const std::string json_string = adapter.serialize(dimensions);
  EXPECT_EQ(json_string, state)
      << "The state mean dimension produced by the "
         "FilterSystemDimensions.deserialize method is incorrect.";
}
/**
 * @brief Test that the KcaStatesJsonAdapter serialize method
 * returns the correct result.
 */
TEST(AdaptersTest, KcaStatesJsonAdapterSerializeTest) {
  // Mock data series to initialise a KCA states object.
  std::vector<std::vector<double>> transition_matrix{
      {1.0011961162353782, 1.0, 0.5}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}
  };
  std::vector<std::vector<double>> transition_covariance{
      {0.12695229227341848, 0, 0}, {0, 0.001, 0}, {0, 0, 0.001}
  };
  std::vector<std::vector<double>> current_state_covariance{
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}
  };
  std::vector<double> current_state_mean{10.288741828687053, 0.0, 0.0};
  std::vector<std::vector<double>> observation_matrix{{1.0, 0.0, 0.0}};
  const double observation_offset{0.0};
  std::vector<double> predicted_state_mean{10.301048359829961, 0.0, 0.0};
  std::vector<std::vector<double>> predicted_state_covariance{
      {0.12695229227341848, 0, 0}, {0, 0.001, 0}, {0, 0, 0.001}
  };

  // Create the KCA states object and set the initial state.
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  kca_states.setTransitionMatrix(transition_matrix);
  kca_states.setTransitionCovariance(transition_covariance);
  kca_states.setCurrentStateCovariance(current_state_covariance);
  kca_states.setCurrentStateMean(current_state_mean);
  kca_states.setObservationMatrix(observation_matrix);
  kca_states.setObservationOffset(observation_offset);
  kca_states.setPredictedStateMean(predicted_state_mean);
  kca_states.setPredictedStateCovariance(predicted_state_covariance);

  // Serialize the KCA states object and create a JSON string.
  const KcaStatesJsonAdapter adapter;
  const std::string state = adapter.serialize(kca_states);

  EXPECT_EQ(
      state, "{\"current_state_covariance\":[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,"
             "0.0,0.0]],\"current_state_mean\":[10.288741828687053,0.0,0.0],"
             "\"observation_matrix\":[[1.0,0.0,0.0]],\"observation_offset\":0."
             "0,\"transition_covariance\":[[0.12695229227341848,0.0,0.0],[0.0,"
             "0.001,0.0],[0.0,0.0,0.001]],\"transition_matrix\":[[1."
             "0011961162353782,1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0]]}"
  ) << "The JSON string produced by the KcaStatesJsonAdapter.serialize "
       "method is incorrect.";
}
/**
 * @brief Test that the KcaStatesJsonAdapter deserialize method
 * returns the correct result.
 */
TEST(AdaptersTest, KcaStatesJsonAdapterDeserializeTest) {
  // Mock data to initialise state objects.
  std::string state =
      "{\"current_state_covariance\":[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,"
      "0.0,0.0]],\"current_state_mean\":[10.288741828687053,0.0,0.0],"
      "\"observation_matrix\":[[1.0,0.0,0.0]],\"observation_offset\":0."
      "0,\"transition_covariance\":[[0.12695229227341848,0.0,0.0],[0.0,"
      "0.001,0.0],[0.0,0.0,0.001]],\"transition_matrix\":[[1."
      "0011961162353782,1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0]]}";
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);

  const KcaStatesJsonAdapter adapter;
  const KcaStates kca_states = adapter.deserialize(state, dimensions);

  // Extract and check the transition, current, observation states from the
  // KCA states object.

  // Transition matrix.
  const std::vector<std::vector<double>> transition_matrix =
      copy_matrix_elements_to_vector(kca_states.getTransitionMatrix());
  const std::vector<std::vector<double>> expected_transition_matrix{
      {1.0011961162353782, 1.0, 0.5}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}
  };
  EXPECT_EQ(transition_matrix, expected_transition_matrix)
      << "The transition matrix was set with invalid or inconsistent values.";

  // Transition covariance.
  const std::vector<std::vector<double>> transition_covariance =
      copy_matrix_elements_to_vector(kca_states.getTransitionCovariance());
  const std::vector<std::vector<double>> expected_transition_covariance{
      {0.12695229227341848, 0, 0}, {0, 0.001, 0}, {0, 0, 0.001}
  };
  EXPECT_EQ(transition_covariance, expected_transition_covariance)
      << "The transition covariance was set with invalid or inconsistent "
         "values.";

  // Current state covariance.
  const std::vector<std::vector<double>> current_state_covariance =
      copy_matrix_elements_to_vector(kca_states.getCurrentStateCovariance());
  const std::vector<std::vector<double>> expected_current_state_covariance{
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}
  };
  EXPECT_EQ(current_state_covariance, expected_current_state_covariance)
      << "The current state covariance matrix was set with invalid or "
         "inconsistent values.";

  // Current state mean.
  std::vector<double> current_state_mean_vector(3);
  const vector<double> current_state_mean = kca_states.getCurrentStateMean();
  std::copy(
      current_state_mean.begin(), current_state_mean.end(),
      current_state_mean_vector.begin()
  );
  const std::vector<double> expected_current_state_mean{
      10.288741828687053, 0.0, 0.0
  };
  EXPECT_EQ(current_state_mean_vector, expected_current_state_mean)
      << "The current state mean vector was set with invalid or inconsistent "
         "values.";

  // Observation matrix.
  const std::vector<std::vector<double>> observation_state =
      copy_matrix_elements_to_vector(kca_states.getObservationMatrix());
  const std::vector<std::vector<double>> expected_observation_matrix{
      {1.0, 0.0, 0.0}
  };
  EXPECT_EQ(observation_state, expected_observation_matrix)
      << "The observation matrix was set with invalid or inconsistent "
         "values.";

  // Observation offset.
  EXPECT_EQ(kca_states.getObservationOffset(), 0.0)
      << "The observation offset value was set to an invalid or inconsistent "
         "value.";
}

/**
 * @test Unparseable JSON previously escaped as a raw
 * nlohmann::json::parse_error because the parse call sat outside the try.
 */
TEST(AdaptersValidationTest, DimensionsDeserializeRejectsUnparseableJson) {
  const FilterSystemDimensionsJsonAdapter adapter;
  EXPECT_THROW(adapter.deserialize("{not json"), json_parse_error)
      << "Unparseable dimensions JSON did not raise json_parse_error.";
}

/**
 * @test Non-integer numeric dimensions previously reached get_to(int&) whose
 * unchecked float->int cast is UB for out-of-range values such as 1e300.
 */
TEST(AdaptersValidationTest, DimensionsDeserializeRejectsNonIntegerFields) {
  const FilterSystemDimensionsJsonAdapter adapter;
  const std::string fractional = R"({"state_mean_dimension": 3.5,
      "state_covariance_rows": 3, "state_covariance_columns": 3,
      "observation_matrix_rows": 1, "observation_matrix_columns": 3,
      "observation_covariance_rows": 1, "observation_covariance_columns": 1,
      "observation_offset": 0.0})";
  EXPECT_THROW(adapter.deserialize(fractional), json_parse_error)
      << "A fractional dimension did not raise json_parse_error.";

  const std::string huge_float = R"({"state_mean_dimension": 1e300,
      "state_covariance_rows": 3, "state_covariance_columns": 3,
      "observation_matrix_rows": 1, "observation_matrix_columns": 3,
      "observation_covariance_rows": 1, "observation_covariance_columns": 1,
      "observation_offset": 0.0})";
  EXPECT_THROW(adapter.deserialize(huge_float), json_parse_error)
      << "A 1e300 dimension did not raise json_parse_error.";

  const std::string stringly = R"({"state_mean_dimension": "3",
      "state_covariance_rows": 3, "state_covariance_columns": 3,
      "observation_matrix_rows": 1, "observation_matrix_columns": 3,
      "observation_covariance_rows": 1, "observation_covariance_columns": 1,
      "observation_offset": 0.0})";
  EXPECT_THROW(adapter.deserialize(stringly), json_parse_error)
      << "A string dimension did not raise json_parse_error.";
}

/**
 * @test Zero, negative, and absurdly large dimensions must be rejected:
 * negatives previously wrapped through size_t matrix constructors into
 * enormous allocations.
 */
TEST(AdaptersValidationTest, DimensionsDeserializeRejectsOutOfRangeFields) {
  const FilterSystemDimensionsJsonAdapter adapter;
  const std::string zero = R"({"state_mean_dimension": 0,
      "state_covariance_rows": 3, "state_covariance_columns": 3,
      "observation_matrix_rows": 1, "observation_matrix_columns": 3,
      "observation_covariance_rows": 1, "observation_covariance_columns": 1,
      "observation_offset": 0.0})";
  EXPECT_THROW(adapter.deserialize(zero), json_parse_error)
      << "A zero dimension did not raise json_parse_error.";

  const std::string negative = R"({"state_mean_dimension": 3,
      "state_covariance_rows": -3, "state_covariance_columns": 3,
      "observation_matrix_rows": 1, "observation_matrix_columns": 3,
      "observation_covariance_rows": 1, "observation_covariance_columns": 1,
      "observation_offset": 0.0})";
  EXPECT_THROW(adapter.deserialize(negative), json_parse_error)
      << "A negative dimension did not raise json_parse_error.";

  const std::string oversized = R"({"state_mean_dimension": 3,
      "state_covariance_rows": 3, "state_covariance_columns": 1025,
      "observation_matrix_rows": 1, "observation_matrix_columns": 3,
      "observation_covariance_rows": 1, "observation_covariance_columns": 1,
      "observation_offset": 0.0})";
  EXPECT_THROW(adapter.deserialize(oversized), json_parse_error)
      << "An out-of-bound dimension did not raise json_parse_error.";
}

/**
 * @test A non-numeric observation_offset must be rejected.
 */
TEST(AdaptersValidationTest, DimensionsDeserializeRejectsNonNumericOffset) {
  const FilterSystemDimensionsJsonAdapter adapter;
  const std::string bad_offset = R"({"state_mean_dimension": 3,
      "state_covariance_rows": 3, "state_covariance_columns": 3,
      "observation_matrix_rows": 1, "observation_matrix_columns": 3,
      "observation_covariance_rows": 1, "observation_covariance_columns": 1,
      "observation_offset": "zero"})";
  EXPECT_THROW(adapter.deserialize(bad_offset), json_parse_error)
      << "A string observation_offset did not raise json_parse_error.";
}

namespace {
// Baseline-valid 3-dim KCA state JSON; each hostile-state test perturbs one
// field. Matches the fixed kinematic scheme dims (3, 3, 3, 1, 3, 1, 1).
const FilterSystemDimensions kKcaSchemeDimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
} // namespace

/**
 * @test Unparseable state JSON previously escaped as a raw
 * nlohmann::json::parse_error because the parse call sat outside the try.
 */
TEST(AdaptersValidationTest, StateDeserializeRejectsUnparseableJson) {
  const KcaStatesJsonAdapter adapter;
  EXPECT_THROW(
      adapter.deserialize("{not json", kKcaSchemeDimensions), json_parse_error
  ) << "Unparseable state JSON did not raise json_parse_error.";
}

/**
 * @test A transition matrix with an extra row previously wrote out of bounds:
 * the move helper trusted the source row count against a 3x3 target.
 */
TEST(AdaptersValidationTest, StateDeserializeRejectsOversizedMatrix) {
  const KcaStatesJsonAdapter adapter;
  const std::string oversized = R"({"current_state_covariance":
      [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]],
      "current_state_mean":[10.28,0.0,0.0],
      "observation_matrix":[[1.0,0.0,0.0]],"observation_offset":0.0,
      "transition_covariance":[[0.1,0.0,0.0],[0.0,0.001,0.0],[0.0,0.0,0.001]],
      "transition_matrix":[[1.0,1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0],
      [9.0,9.0,9.0]]})";
  EXPECT_THROW(
      adapter.deserialize(oversized, kKcaSchemeDimensions), json_parse_error
  ) << "An oversized transition_matrix did not raise json_parse_error.";
}

/**
 * @test A ragged inner row must be rejected before any write.
 */
TEST(AdaptersValidationTest, StateDeserializeRejectsRaggedMatrix) {
  const KcaStatesJsonAdapter adapter;
  const std::string ragged = R"({"current_state_covariance":
      [[0.0,0.0,0.0],[0.0,0.0,0.0,7.0],[0.0,0.0,0.0]],
      "current_state_mean":[10.28,0.0,0.0],
      "observation_matrix":[[1.0,0.0,0.0]],"observation_offset":0.0,
      "transition_covariance":[[0.1,0.0,0.0],[0.0,0.001,0.0],[0.0,0.0,0.001]],
      "transition_matrix":[[1.0,1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0]]})";
  EXPECT_THROW(
      adapter.deserialize(ragged, kKcaSchemeDimensions), json_parse_error
  ) << "A ragged current_state_covariance did not raise json_parse_error.";
}

/**
 * @test An undersized state mean must be rejected (uninitialised remainder
 * otherwise).
 */
TEST(AdaptersValidationTest, StateDeserializeRejectsWrongLengthMean) {
  const KcaStatesJsonAdapter adapter;
  const std::string short_mean = R"({"current_state_covariance":
      [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]],
      "current_state_mean":[10.28,0.0],
      "observation_matrix":[[1.0,0.0,0.0]],"observation_offset":0.0,
      "transition_covariance":[[0.1,0.0,0.0],[0.0,0.001,0.0],[0.0,0.0,0.001]],
      "transition_matrix":[[1.0,1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0]]})";
  EXPECT_THROW(
      adapter.deserialize(short_mean, kKcaSchemeDimensions), json_parse_error
  ) << "A 2-element current_state_mean did not raise json_parse_error.";
}

/**
 * @test An observation matrix whose shape disagrees with the dimensions input
 * must be rejected.
 */
TEST(AdaptersValidationTest, StateDeserializeRejectsWrongShapeObservationMatrix) {
  const KcaStatesJsonAdapter adapter;
  const std::string wrong_shape = R"({"current_state_covariance":
      [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]],
      "current_state_mean":[10.28,0.0,0.0],
      "observation_matrix":[[1.0,0.0,0.0],[2.0,0.0,0.0]],
      "observation_offset":0.0,
      "transition_covariance":[[0.1,0.0,0.0],[0.0,0.001,0.0],[0.0,0.0,0.001]],
      "transition_matrix":[[1.0,1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0]]})";
  EXPECT_THROW(
      adapter.deserialize(wrong_shape, kKcaSchemeDimensions), json_parse_error
  ) << "A 2x3 observation_matrix did not raise json_parse_error for 1x3 "
       "dimensions.";
}

/**
 * @test A non-numeric observation_offset in the state JSON must be rejected.
 */
TEST(AdaptersValidationTest, StateDeserializeRejectsNonNumericOffset) {
  const KcaStatesJsonAdapter adapter;
  const std::string bad_offset = R"({"current_state_covariance":
      [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]],
      "current_state_mean":[10.28,0.0,0.0],
      "observation_matrix":[[1.0,0.0,0.0]],"observation_offset":"zero",
      "transition_covariance":[[0.1,0.0,0.0],[0.0,0.001,0.0],[0.0,0.0,0.001]],
      "transition_matrix":[[1.0,1.0,0.5],[0.0,1.0,1.0],[0.0,0.0,1.0]]})";
  EXPECT_THROW(
      adapter.deserialize(bad_offset, kKcaSchemeDimensions), json_parse_error
  ) << "A string observation_offset did not raise json_parse_error.";
}
