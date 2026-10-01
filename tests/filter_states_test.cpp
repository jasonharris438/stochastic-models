#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/kalman_filter/states.h"
#include "stochastic_models/kalman_filter/states_exceptions.h"
#include "stochastic_models/kalman_filter/type_conversion.h"
#include "stochastic_models/numeric_utils/linalg.h"
#include "support/assertions.h"
#include "support/expected_values.h"

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <stdexcept>
#include <type_traits>
#include <vector>

using namespace boost::numeric::ublas;

constexpr double kalman_tolerance = 1e-12;

using test_support::expectMatrixNear;
using test_support::expectVectorNear;
/**
 * @brief Function to create a transition matrix for testing.
 * @return std::vector<std::vector<double>>. The transition matrix.
 */
std::vector<std::vector<double>> transition_matrix() {
  return {{1.000295, 1.0, 0.5}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}};
}
/**
 * @brief Function to create a transition matrix in a boost object for testing.
 *
 * Creates a matrix from a std::vector and then moves the values into the boost
 * matrix.
 *
 * @return matrix<double>. The boost matrix.
 */
const matrix<double> create_boost_transition_matrix() {
  // Create a std::vector first.
  std::vector<std::vector<double>> std_vectors = transition_matrix();

  // Allocate a boost matrix and then move the values from the std::vector.
  matrix<double> boost_matrix(std_vectors.size(), std_vectors[0].size());
  for (int i{0}; i < std_vectors.size(); i++) {
    std::move(
        std_vectors.at(i).begin(), std_vectors.at(i).end(),
        row(boost_matrix, i).begin()
    );
  }
  return boost_matrix;
}
/**
 * @brief Function to create a current covariance matrix for testing as a
 * std::vector.
 * @return std::vector<std::vector<double>>. The current covariance matrix as a
 * std::vector.
 */
std::vector<std::vector<double>> current_covariance() {
  return {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
}
/**
 * @brief Function to create a current covariance matrix in a boost object for
 * testing.
 *
 * Creates a matrix from a std::vector and then moves the values into the boost
 * matrix.
 *
 * @return matrix<double>. The boost matrix.
 */
const matrix<double> create_boost_current_covariance_matrix() {
  std::vector<std::vector<double>> std_vectors = current_covariance();
  matrix<double> current_cov(std_vectors.size(), std_vectors[0].size());
  for (int i{0}; i < std_vectors.size(); i++) {
    std::move(
        std_vectors.at(i).begin(), std_vectors.at(i).end(),
        row(current_cov, i).begin()
    );
  }
  return current_cov;
}
/**
 * @brief Function to create a boost current mean vector for testing.
 *
 * Creates a std::vector and then moves the values into the boost vector.
 *
 * @return std::vector<double>. A boost vector.
 */
const vector<double> create_current_mean_vector() {
  std::vector<double> std_vector{1.330593, 0.0, 0.0};
  vector<double> boost_vector(std_vector.size());
  std::move(std_vector.begin(), std_vector.end(), boost_vector.begin());
  return boost_vector;
}
/**
 * @brief Test that the PredictedState.calculateCovariance method returns
 * the correct result.
 */
TEST(KalmanFilterTest, PredictedStateCalculateCovarianceTest) {
  const matrix<double> boost_matrix = create_boost_transition_matrix();
  const PredictedState matrix_test(boost_matrix);

  std::vector<std::vector<double>> transition_covariance_std_vector{
      {0.013744, 0.0, 0.0}, {0.0, 0.001, 0.0}, {0.0, 0.0, 0.001}
  };
  const matrix<double> transition_covariance =
      create_boost_matrix_from_vectors(transition_covariance_std_vector);

  const matrix<double> current_cov = create_boost_current_covariance_matrix();

  const matrix<double> result =
      matrix_test.calculateCovariance(current_cov, transition_covariance);
  const std::vector<std::vector<double>> result_vector =
      copy_matrix_elements_to_vector(result);

  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      result_vector,
      expected::kalman_filter_test::
          predicted_state_calculate_covariance_test_predicted_state_covariance,
      kalman_tolerance, "The covariance matrix is not calculated correctly."
  ));
}
/**
 * @brief Test that the PredictedState.calculateMean method returns the correct
 * result.
 */
TEST(KalmanFilterTest, PredictedStateCalculateMeanTest) {
  const matrix<double> boost_matrix = create_boost_transition_matrix();
  const PredictedState matrix_test(boost_matrix);

  const vector<double> boost_vector = create_current_mean_vector();
  const vector<double> result = matrix_test.calculateMean(boost_vector);
  const std::vector<double> result_vector(result.begin(), result.end());

  ASSERT_NO_FATAL_FAILURE(expectVectorNear(
      result_vector,
      expected::kalman_filter_test::
          predicted_state_calculate_mean_test_predicted_state_mean,
      kalman_tolerance, "The mean vector is not calculated correctly."
  ));
}
/**
 * @brief Test that the PredictedObservation.calculateMean method returns the
 * correct result.
 */
TEST(KalmanFilterTest, PredictedObservationCalculateMeanTest) {
  matrix<double> observation_matrix = matrix<double>(1, 3);
  observation_matrix(0, 0) = 1;
  observation_matrix(0, 1) = 0;
  observation_matrix(0, 2) = 0;

  const double observation_offset = 0;
  const PredictedObservation predicted_observation(
      observation_matrix, observation_offset
  );

  vector<double> predicted_state = vector<double>(3);
  predicted_state(0) = 1.330986;
  predicted_state(1) = 0;
  predicted_state(2) = 0;

  const vector<double> result =
      predicted_observation.calculateMean(predicted_state);
  const std::vector<double> result_vector(result.begin(), result.end());

  ASSERT_NO_FATAL_FAILURE(expectVectorNear(
      result_vector,
      expected::kalman_filter_test::
          predicted_observation_calculate_mean_test_predicted_observation_mean,
      kalman_tolerance, "The mean vector is not calculated correctly."
  ));
}
/**
 * @brief Test that the PredictedObservation.calculateCovariance method returns
 * the correct result.
 */
TEST(KalmanFilterTest, PredictedObservationCalculateCovarianceTest) {
  // Create example observation matrix.
  matrix<double> observation_matrix = matrix<double>(1, 3);
  observation_matrix(0, 0) = 1;
  observation_matrix(0, 1) = 0;
  observation_matrix(0, 2) = 0;

  // Example observation offset.
  const double observation_offset = 0;
  const PredictedObservation predicted_observation(
      observation_matrix, observation_offset
  );

  // Predicted state covariance from a std::vector.
  const std::vector<std::vector<double>> covariance_std_vector{
      {0.013744, 0.0, 0.0}, {0.0, 0.001, 0.0}, {0.0, 0.0, 0.001}
  };
  const matrix<double> predicted_state_covariance =
      create_boost_matrix_from_vectors(covariance_std_vector);

  // Sigma parameter from SDE.
  const double sigma_param{0.00687526};

  // Generate predicted observation covariance from variables created.
  const matrix<double> result = predicted_observation.calculateCovariance(
      predicted_state_covariance, sigma_param
  );
  const std::vector<std::vector<double>> result_vec =
      copy_matrix_elements_to_vector(result);

  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      result_vec,
      expected::kalman_filter_test::
          predicted_observation_calculate_covariance_test_predicted_observation_covariance,
      kalman_tolerance,
      "The predicted observation covariance is not calculated correctly."
  ));
}
/**
 * @brief Test that the PredictedObservation.calculateKalmanGain method returns
 * the correct result.
 */
TEST(KalmanFilterTest, PredictedObservationCalculateKalmanGainTest) {
  // Create example observation matrix.
  matrix<double> observation_matrix = matrix<double>(1, 3);
  observation_matrix(0, 0) = 1;
  observation_matrix(0, 1) = 0;
  observation_matrix(0, 2) = 0;

  // Example observation offset.
  const double observation_offset = 0;
  const PredictedObservation predicted_observation(
      observation_matrix, observation_offset
  );

  // Predicted observation covariance.
  matrix<double> predicted_observation_covariance = matrix<double>(1, 1);
  predicted_observation_covariance(0, 0) = 0.0137917;

  // Predicted state covariance.
  const std::vector<std::vector<double>> covariance_vector{
      {0.013744, 0.0, 0.0}, {0.0, 0.001, 0.0}, {0.0, 0.0, 0.001}
  };
  const matrix<double> predicted_state_covariance =
      create_boost_matrix_from_vectors(covariance_vector);

  // Generate the Kalman gain from variables created.
  const BoostMatrixInverter matrix_inverter;
  const matrix<double> result = predicted_observation.calculateKalmanGain(
      predicted_state_covariance, predicted_observation_covariance,
      matrix_inverter
  );
  const std::vector<std::vector<double>> result_vec =
      copy_matrix_elements_to_vector(result);

  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      result_vec,
      expected::kalman_filter_test::
          predicted_observation_calculate_kalman_gain_test_kalman_gain,
      kalman_tolerance, "The Kalman gain is not calculated correctly."
  ));
}
/**
 * @brief Test that the CurrentState.calculateMean method returns the
 * correct result.
 */
TEST(KalmanFilterTest, CurrentStateCalculateMeanTest) {
  const double innovation{-0.02018567};

  vector<double> predicted_state_mean = vector<double>(3);
  predicted_state_mean(0) = 1.330986;
  predicted_state_mean(1) = 0;
  predicted_state_mean(2) = 0;

  vector<double> kalman_gain(3);
  kalman_gain(0) = 0.99657263;
  kalman_gain(1) = 0;
  kalman_gain(2) = 0;

  const CurrentState current_state;

  const vector<double> result = current_state.calculateMean(
      predicted_state_mean, kalman_gain, innovation
  );
  const std::vector<double> result_vector(result.begin(), result.end());

  ASSERT_NO_FATAL_FAILURE(expectVectorNear(
      result_vector,
      expected::kalman_filter_test::
          current_state_calculate_mean_test_current_state_mean,
      kalman_tolerance, "The mean vector is not calculated correctly."
  ));
}

/**
 * @brief Test that the CurrentState.calculateCovariance method returns the
 * correct result.
 */
TEST(KalmanFilterTest, CurrentStateCalculateCovarianceTest) {
  // Create example observation matrix.
  matrix<double> observation_matrix = matrix<double>(1, 3);
  observation_matrix(0, 0) = 1;
  observation_matrix(0, 1) = 0;
  observation_matrix(0, 2) = 0;

  // Create example kalman gain.
  matrix<double> kalman_gain(3, 1);
  kalman_gain(0, 0) = 0.99657263;
  kalman_gain(1, 0) = 0;
  kalman_gain(2, 0) = 0;

  // Predicted state covariance.
  const std::vector<std::vector<double>> covariance_vector{
      {0.013744, 0.0, 0.0}, {0.0, 0.001, 0.0}, {0.0, 0.0, 0.001}
  };
  const matrix<double> predicted_state_covariance =
      create_boost_matrix_from_vectors(covariance_vector);

  const CurrentState current_state;

  const matrix<double> result = current_state.calculateCovariance(
      predicted_state_covariance, observation_matrix, kalman_gain
  );

  const std::vector<std::vector<double>> result_vec =
      copy_matrix_elements_to_vector(result);

  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      result_vec,
      expected::kalman_filter_test::
          current_state_calculate_covariance_test_current_state_covariance,
      kalman_tolerance, "The covariance matrix is not calculated correctly."
  ));
}
/**
 * @brief Test that the KcaStates setInitialState sets a KcaStates instance to
 * the correct initial state.
 */
TEST(KalmanFilterStateTest, KcaStatessetInitialStateTest) {
  // This test checks all components that are set by the setInitialState
  // method of the KcaStates class.

  // Mock data series to initialise the KCA states object.
  const std::vector<double> data_series{10.51255, 10.51985, 10.52405, 10.4656,
                                        10.47,    10.5403,  10.4425,  10.3087,
                                        10.1994,  10.1839,  10.24645, 10.1795,
                                        10.21715, 10.14995, 10.194,   10.22505,
                                        10.27325, 10.25095, 10.30575, 10.27645};
  const double h{1.0};
  const double q{0.001};

  // Create the KCA states object and check the initial state is set
  // correctly.
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  kca_states.setInitialState(data_series, h, q);
  EXPECT_TRUE(kca_states.isInitialised())
      << "The KCA states object must correctly indicate whether it is "
         "initialised.";

  // Extract and check the transition, current, observation states from the
  // KCA states object.

  // Transition matrix.
  const std::vector<std::vector<double>> transition_matrix =
      copy_matrix_elements_to_vector(kca_states.getTransitionMatrix());
  const std::vector<std::vector<double>> expected_transition_matrix{
      {1.0, 1.0, 0.5}, {0.0, 1.0, 1.0}, {0.0, 0.0, 1.0}
  };
  EXPECT_EQ(transition_matrix, expected_transition_matrix)
      << "The transition matrix was set with invalid or inconsistent values.";

  // Transition covariance.
  const std::vector<std::vector<double>> transition_covariance =
      copy_matrix_elements_to_vector(kca_states.getTransitionCovariance());
  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      transition_covariance,
      expected::kalman_filter_state_test::
          kca_statesset_initial_state_test_transition_covariance,
      kalman_tolerance,
      "The transition covariance was set with invalid or inconsistent values."
  ));

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
  const vector<double> current_state_mean = kca_states.getCurrentStateMean();
  const std::vector<double> current_state_mean_vector(
      current_state_mean.begin(), current_state_mean.end()
  );
  const std::vector<double> expected_current_state_mean{
      data_series.back(), 0.0, 0.0
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
 * @brief Test that the KcaStates updatePredictedState sets a KcaStates instance
 * to the correct state with valid prior values.
 */
TEST(KalmanFilterStateTest, KcaStatesupdatePredictedStateTest) {
  // This test checks all components that are updated by the
  // updatePredictedState method of the KcaStates class.

  // Mock data series to initialise the KCA states object.
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

  // Create the KCA states object and check the initial state is set
  // correctly.
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  kca_states.setTransitionMatrix(transition_matrix);
  kca_states.setTransitionCovariance(transition_covariance);
  kca_states.setCurrentStateCovariance(current_state_covariance);
  kca_states.setCurrentStateMean(current_state_mean);
  kca_states.setObservationMatrix(observation_matrix);
  kca_states.setObservationOffset(observation_offset);
  kca_states.setInitialized();

  // Extract and check the transition, current, observation states from the
  // KCA states object.
  kca_states.updatePredictedState();
  EXPECT_TRUE(kca_states.arePriorsValid())
      << "The KCA states object must correctly indicate whether the prior "
         "state is valid.";

  // Predicted state mean.
  const vector<double> predicted_state_mean =
      kca_states.getPredictedStateMean();
  const std::vector<double> predicted_state_mean_vector(
      predicted_state_mean.begin(), predicted_state_mean.end()
  );
  ASSERT_NO_FATAL_FAILURE(expectVectorNear(
      predicted_state_mean_vector,
      expected::kalman_filter_state_test::
          kca_statesupdate_predicted_state_test_predicted_state_mean,
      kalman_tolerance,
      "The predicted state mean vector was set with invalid or "
      "inconsistent values."
  ));

  // Predicted state covariance.
  const std::vector<std::vector<double>> predicted_state_covariance =
      copy_matrix_elements_to_vector(kca_states.getPredictedStateCovariance());
  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      predicted_state_covariance,
      expected::kalman_filter_state_test::
          kca_statesupdate_predicted_state_test_predicted_state_covariance,
      kalman_tolerance,
      "The predicted state covariance was set with invalid or "
      "inconsistent values."
  ));
}
/**
 * @brief Test that the KcaStates updateCurrentState sets a KcaStates instance
 * to the correct state with valid posterior values.
 */
TEST(KalmanFilterStateTest, KcaStatesupdateCurrentStateTest) {
  // This test checks all components that are updated by the
  // updateCurrentState method of the KcaStates class.

  // Mock data series to initialise the KCA states object.
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
  std::vector<double> predicted_observation_mean{10.301};

  // Create the KCA states object and check the initial state is set
  // correctly.
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
  kca_states.setPredictedObservationMean(predicted_observation_mean);
  kca_states.setInitialized();
  kca_states.setPriorsTrue();

  // Extract and check the predicted observation and current states from the
  // KCA states object.
  kca_states.updateCurrentState(10.3, 0.1);
  EXPECT_FALSE(kca_states.arePriorsValid())
      << "The KCA states object must correctly invalidate the prior state "
         "when updateCurrentState is called.";

  // Predicted observation mean.
  const vector<double> updated_predicted_observation_mean =
      kca_states.getPredictedObservationMean();
  const std::vector<double> predicted_observation_mean_vector(
      updated_predicted_observation_mean.begin(),
      updated_predicted_observation_mean.end()
  );
  ASSERT_NO_FATAL_FAILURE(expectVectorNear(
      predicted_observation_mean_vector,
      expected::kalman_filter_state_test::
          kca_statesupdate_current_state_test_predicted_observation_mean,
      kalman_tolerance,
      "The predicted observation mean vector was set with invalid or "
      "inconsistent values."
  ));

  // Predicted observation covariance.
  const std::vector<std::vector<double>>
      predicted_observation_covariance_vector = copy_matrix_elements_to_vector(
          kca_states.getPredictedObservationCovariance()
      );
  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      predicted_observation_covariance_vector,
      expected::kalman_filter_state_test::
          kca_statesupdate_current_state_test_predicted_observation_covariance,
      kalman_tolerance,
      "The predicted observation covariance was set with invalid or "
      "inconsistent values."
  ));

  // Current state mean.
  const vector<double> new_current_state_mean =
      kca_states.getCurrentStateMean();
  const std::vector<double> new_current_state_mean_vector(
      new_current_state_mean.begin(), new_current_state_mean.end()
  );
  ASSERT_NO_FATAL_FAILURE(expectVectorNear(
      new_current_state_mean_vector,
      expected::kalman_filter_state_test::
          kca_statesupdate_current_state_test_current_state_mean,
      kalman_tolerance,
      "The current state mean vector was set with invalid or "
      "inconsistent values."
  ));

  // Current state covariance.
  const std::vector<std::vector<double>> new_current_state_covariance =
      copy_matrix_elements_to_vector(kca_states.getCurrentStateCovariance());
  ASSERT_NO_FATAL_FAILURE(expectMatrixNear(
      new_current_state_covariance,
      expected::kalman_filter_state_test::
          kca_statesupdate_current_state_test_current_state_covariance,
      kalman_tolerance,
      "The current state covariance was set with invalid or "
      "inconsistent values."
  ));
}

/**
 * @test A source matrix with more rows than the target previously wrote out
 * of bounds through unchecked uBLAS row proxies. The move helper must verify
 * exact row/column counts.
 */
TEST(FilterStatesValidationTest, setTransitionMatrixRejectsRowCountMismatch) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  std::vector<std::vector<double>> oversized{
      {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {9.0, 9.0, 9.0}
  };
  EXPECT_THROW(kca_states.setTransitionMatrix(oversized), filter_shape_mismatch)
      << "setTransitionMatrix accepted a 4-row source for a 3x3 target.";
}

/**
 * @test A ragged source row longer than the target's column count previously
 * overran the row. It must be rejected.
 */
TEST(FilterStatesValidationTest, setTransitionMatrixRejectsRaggedRows) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  std::vector<std::vector<double>> ragged{
      {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 5.0}, {0.0, 0.0, 1.0}
  };
  EXPECT_THROW(kca_states.setTransitionMatrix(ragged), filter_shape_mismatch)
      << "setTransitionMatrix accepted a ragged source row.";
}

/**
 * @test An oversized source vector must be rejected by the vector move helper.
 */
TEST(FilterStatesValidationTest, setCurrentStateMeanRejectsLengthMismatch) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  std::vector<double> oversized{1.0, 2.0, 3.0, 4.0};
  EXPECT_THROW(kca_states.setCurrentStateMean(oversized), filter_shape_mismatch)
      << "setCurrentStateMean accepted a 4-element source for a 3-element "
         "target.";
}

/**
 * @test The uBLAS overload of setCurrentStateMean must reject a source that is
 * longer than the target instead of writing past it.
 */
TEST(
    FilterStatesValidationTest, setCurrentStateMeanUblasRejectsLengthMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const vector<double> oversized = scalar_vector<double>(4, 1.0);
  EXPECT_THROW(kca_states.setCurrentStateMean(oversized), filter_shape_mismatch)
      << "setCurrentStateMean accepted a 4-element uBLAS source for a "
         "3-element target.";
}

/**
 * @test The uBLAS overload of setPredictedStateMean must reject a source that
 * is longer than the target.
 */
TEST(
    FilterStatesValidationTest, setPredictedStateMeanUblasRejectsLengthMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const vector<double> oversized = scalar_vector<double>(4, 1.0);
  EXPECT_THROW(
      kca_states.setPredictedStateMean(oversized), filter_shape_mismatch
  ) << "setPredictedStateMean accepted a 4-element uBLAS source for a "
       "3-element target.";
}

/**
 * @test A shorter source must also be rejected. Previously a length-2 source
 * into a length-1 target wrote one element past the heap allocation.
 */
TEST(
    FilterStatesValidationTest,
    setPredictedObservationMeanUblasRejectsLengthMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const vector<double> undersized = scalar_vector<double>(2, 1.0);
  EXPECT_THROW(
      kca_states.setPredictedObservationMean(undersized), filter_shape_mismatch
  ) << "setPredictedObservationMean accepted a 2-element uBLAS source for a "
       "1-element target.";
}

/**
 * @test The dimension and shape exception types must be catchable as
 * std::invalid_argument, so existing catch sites keep working.
 */
TEST(FilterStatesValidationTest, shapeExceptionsDeriveFromInvalidArgument) {
  static_assert(
      std::is_base_of_v<std::invalid_argument, invalid_filter_dimensions>
  );
  static_assert(
      std::is_base_of_v<std::invalid_argument, filter_shape_mismatch>
  );
  const invalid_filter_dimensions dimension_error{"dimension"};
  const filter_shape_mismatch shape_error{"shape"};
  EXPECT_STREQ(dimension_error.what(), "dimension")
      << "invalid_filter_dimensions did not keep its message.";
  EXPECT_STREQ(shape_error.what(), "shape")
      << "filter_shape_mismatch did not keep its message.";
}

/**
 * @test A zero dimension must be rejected at construction.
 */
TEST(FilterStatesValidationTest, dimensionsRejectZeroDimension) {
  EXPECT_THROW(
      FilterSystemDimensions(0, 3, 3, 1, 3, 1, 1, 0.0),
      invalid_filter_dimensions
  ) << "FilterSystemDimensions accepted a zero state mean dimension.";
}

/**
 * @test A dimension above the shared bound must be rejected at construction.
 */
TEST(FilterStatesValidationTest, dimensionsRejectDimensionAboveBound) {
  const std::size_t above_bound = FilterSystemDimensions::max_dimension + 1;
  EXPECT_THROW(
      FilterSystemDimensions(
          above_bound, above_bound, above_bound, 1, above_bound, 1, 1, 0.0
      ),
      invalid_filter_dimensions
  ) << "FilterSystemDimensions accepted a dimension above max_dimension.";
}

/**
 * @test A negative literal converts to a huge std::size_t in paren-init.
 * Previously it wrapped into a one-element matrix that reported both
 * dimensions as 2^64 - 1. It must fail the upper bound.
 */
TEST(FilterStatesValidationTest, dimensionsRejectNegativeDimension) {
  EXPECT_THROW(
      FilterSystemDimensions(3, -1, -1, 1, 3, 1, 1, 0.0),
      invalid_filter_dimensions
  ) << "FilterSystemDimensions accepted a negative state covariance "
       "dimension.";
}

/**
 * @test The KCA scheme constructs and exposes each dimension through an
 * accessor.
 */
TEST(FilterStatesValidationTest, dimensionsAcceptKcaScheme) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.5);
  EXPECT_EQ(dimensions.getStateMeanDimension(), 3u)
      << "getStateMeanDimension returned the wrong dimension.";
  EXPECT_EQ(dimensions.getStateCovarianceRows(), 3u)
      << "getStateCovarianceRows returned the wrong dimension.";
  EXPECT_EQ(dimensions.getStateCovarianceColumns(), 3u)
      << "getStateCovarianceColumns returned the wrong dimension.";
  EXPECT_EQ(dimensions.getObservationMatrixRows(), 1u)
      << "getObservationMatrixRows returned the wrong dimension.";
  EXPECT_EQ(dimensions.getObservationMatrixColumns(), 3u)
      << "getObservationMatrixColumns returned the wrong dimension.";
  EXPECT_EQ(dimensions.getObservationCovarianceRows(), 1u)
      << "getObservationCovarianceRows returned the wrong dimension.";
  EXPECT_EQ(dimensions.getObservationCovarianceColumns(), 1u)
      << "getObservationCovarianceColumns returned the wrong dimension.";
  EXPECT_EQ(dimensions.getObservationOffset(), 0.5)
      << "getObservationOffset returned the wrong value.";
}

/**
 * @test The four state dimensions must be equal. The filter algebra
 * multiplies the transition matrix into the state mean, so a mismatch is
 * unrunnable.
 */
TEST(FilterStatesValidationTest, dimensionsRejectUnequalStateDimensions) {
  EXPECT_THROW(
      FilterSystemDimensions(3, 1, 3, 1, 3, 1, 1, 0.0),
      invalid_filter_dimensions
  ) << "FilterSystemDimensions accepted state covariance rows of 1 with a "
       "state mean dimension of 3.";
}

/**
 * @test The observation covariance must be square in the observation matrix
 * row count.
 */
TEST(FilterStatesValidationTest, dimensionsRejectUnequalObservationDimensions) {
  EXPECT_THROW(
      FilterSystemDimensions(3, 3, 3, 1, 3, 2, 2, 0.0),
      invalid_filter_dimensions
  ) << "FilterSystemDimensions accepted a 2x2 observation covariance with a "
       "1-row observation matrix.";
}

/**
 * @test A consistent set that is not the KCA scheme must construct, so the
 * type stays generic in the state dimension.
 */
TEST(FilterStatesValidationTest, dimensionsAcceptConsistentNonKcaScheme) {
  const FilterSystemDimensions dimensions(4, 4, 4, 1, 4, 1, 1, 1.0);
  EXPECT_EQ(dimensions.getStateMeanDimension(), 4u)
      << "A consistent 4-state dimension set did not construct.";
}

/**
 * @test A consistent vector-observation set must construct. The scalar
 * observation rule belongs to KcaStates, not to the dimensions type.
 */
TEST(FilterStatesValidationTest, dimensionsAcceptVectorObservationScheme) {
  const FilterSystemDimensions dimensions(3, 3, 3, 2, 3, 2, 2, 0.0);
  EXPECT_EQ(dimensions.getObservationMatrixRows(), 2u)
      << "A consistent 2-row observation set did not construct.";
}

/**
 * @test The filter takes one scalar observation per update, so KcaStates must
 * reject a dimension set with more than one observation row. This is the
 * dimension set behind the first reported heap write.
 */
TEST(FilterStatesValidationTest, kcaStatesRejectVectorObservationDimensions) {
  const FilterSystemDimensions dimensions(3, 3, 3, 2, 3, 2, 2, 0.0);
  EXPECT_THROW(KcaStates rejected(dimensions), invalid_filter_dimensions)
      << "KcaStates accepted a 2-row observation matrix.";
}

/**
 * @test An empty series has no last value to seed the state mean. It must be
 * rejected with the library's observation-count exception, not with an STL
 * range error.
 */
TEST(FilterStatesValidationTest, setInitialStateRejectsEmptySeries) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  EXPECT_THROW(
      kca_states.setInitialState({}, 1.0, 0.001), InvalidNumberObservationsError
  ) << "setInitialState accepted an empty data series.";
}

/**
 * @test The uBLAS overload of setTransitionMatrix must reject a source with
 * more rows than the target.
 */
TEST(FilterStatesValidationTest, setTransitionMatrixUblasRejectsRowMismatch) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const matrix<double> oversized = scalar_matrix<double>(4, 3, 1.0);
  EXPECT_THROW(kca_states.setTransitionMatrix(oversized), filter_shape_mismatch)
      << "setTransitionMatrix accepted a 4x3 uBLAS source for a 3x3 target.";
}

/**
 * @test The uBLAS overload of setTransitionMatrix must reject a source with
 * more columns than the target.
 */
TEST(
    FilterStatesValidationTest, setTransitionMatrixUblasRejectsColumnMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const matrix<double> oversized = scalar_matrix<double>(3, 4, 1.0);
  EXPECT_THROW(kca_states.setTransitionMatrix(oversized), filter_shape_mismatch)
      << "setTransitionMatrix accepted a 3x4 uBLAS source for a 3x3 target.";
}

/**
 * @test The uBLAS overload of setTransitionCovariance must reject a row
 * mismatch.
 */
TEST(
    FilterStatesValidationTest, setTransitionCovarianceUblasRejectsRowMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const matrix<double> oversized = scalar_matrix<double>(4, 3, 1.0);
  EXPECT_THROW(
      kca_states.setTransitionCovariance(oversized), filter_shape_mismatch
  ) << "setTransitionCovariance accepted a 4x3 uBLAS source for a 3x3 "
       "target.";
}

/**
 * @test The uBLAS overload of setCurrentStateCovariance must reject a row
 * mismatch.
 */
TEST(
    FilterStatesValidationTest, setCurrentStateCovarianceUblasRejectsRowMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const matrix<double> oversized = scalar_matrix<double>(4, 3, 1.0);
  EXPECT_THROW(
      kca_states.setCurrentStateCovariance(oversized), filter_shape_mismatch
  ) << "setCurrentStateCovariance accepted a 4x3 uBLAS source for a 3x3 "
       "target.";
}

/**
 * @test The uBLAS overload of setPredictedStateCovariance must reject a row
 * mismatch.
 */
TEST(
    FilterStatesValidationTest,
    setPredictedStateCovarianceUblasRejectsRowMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const matrix<double> oversized = scalar_matrix<double>(4, 3, 1.0);
  EXPECT_THROW(
      kca_states.setPredictedStateCovariance(oversized), filter_shape_mismatch
  ) << "setPredictedStateCovariance accepted a 4x3 uBLAS source for a 3x3 "
       "target.";
}

/**
 * @test The uBLAS overload of setObservationMatrix must reject a row mismatch
 * against its 1x3 target.
 */
TEST(FilterStatesValidationTest, setObservationMatrixUblasRejectsRowMismatch) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const matrix<double> oversized = scalar_matrix<double>(2, 3, 1.0);
  EXPECT_THROW(
      kca_states.setObservationMatrix(oversized), filter_shape_mismatch
  ) << "setObservationMatrix accepted a 2x3 uBLAS source for a 1x3 target.";
}

/**
 * @test The uBLAS overload of setPredictedObservationCovariance must reject a
 * shape mismatch against its 1x1 target.
 */
TEST(
    FilterStatesValidationTest,
    setPredictedObservationCovarianceUblasRejectsShapeMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const matrix<double> oversized = scalar_matrix<double>(2, 2, 1.0);
  EXPECT_THROW(
      kca_states.setPredictedObservationCovariance(oversized),
      filter_shape_mismatch
  ) << "setPredictedObservationCovariance accepted a 2x2 uBLAS source for a "
       "1x1 target.";
}

/**
 * @test A matching uBLAS source must be copied element for element and must
 * leave the target shape unchanged.
 */
TEST(FilterStatesValidationTest, setTransitionMatrixUblasCopiesMatchingSource) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const matrix<double> source = create_boost_transition_matrix();
  kca_states.setTransitionMatrix(source);
  const matrix<double>& stored = kca_states.getTransitionMatrix();
  EXPECT_EQ(stored.size1(), 3u)
      << "setTransitionMatrix changed the target row count.";
  EXPECT_EQ(stored.size2(), 3u)
      << "setTransitionMatrix changed the target column count.";
  EXPECT_EQ(
      copy_matrix_elements_to_vector(stored),
      copy_matrix_elements_to_vector(source)
  ) << "setTransitionMatrix did not copy every element of the source.";
}

/**
 * @test KcaStates writes a fixed three-state kinematic scheme, so it must
 * reject a consistent four-state dimension set at construction.
 */
TEST(FilterStatesValidationTest, kcaStatesRejectFourStateDimensions) {
  const FilterSystemDimensions dimensions(4, 4, 4, 1, 4, 1, 1, 0.0);
  EXPECT_THROW(KcaStates rejected(dimensions), invalid_filter_dimensions)
      << "KcaStates accepted a 4-state dimension set.";
}

/**
 * @test A PriorState built from a dimensions object must size every buffer
 * from that object, so no raw size can reach uBLAS.
 */
TEST(FilterStatesValidationTest, priorStateSizesBuffersFromDimensions) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.5);
  const PriorState prior_state(dimensions);
  EXPECT_EQ(prior_state.predicted_observation_mean.size(), 1u)
      << "predicted_observation_mean was not sized by the observation rows.";
  EXPECT_EQ(prior_state.predicted_state_mean.size(), 3u)
      << "predicted_state_mean was not sized by the state mean dimension.";
  EXPECT_EQ(prior_state.predicted_observation_covariance.size1(), 1u)
      << "predicted_observation_covariance rows were not sized by the "
         "dimensions.";
  EXPECT_EQ(prior_state.predicted_state_covariance.size2(), 3u)
      << "predicted_state_covariance columns were not sized by the dimensions.";
  EXPECT_EQ(prior_state.observation_matrix.size2(), 3u)
      << "observation_matrix columns were not sized by the dimensions.";
  EXPECT_DOUBLE_EQ(prior_state.observation_offset, 0.5)
      << "observation_offset was not copied from the dimensions.";
}

/**
 * @test A PosteriorState built from a dimensions object must size both
 * buffers from that object.
 */
TEST(FilterStatesValidationTest, posteriorStateSizesBuffersFromDimensions) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  const PosteriorState posterior_state(dimensions);
  EXPECT_EQ(posterior_state.current_state_mean.size(), 3u)
      << "current_state_mean was not sized by the state mean dimension.";
  EXPECT_EQ(posterior_state.current_state_covariance.size1(), 3u)
      << "current_state_covariance rows were not sized by the dimensions.";
  EXPECT_EQ(posterior_state.current_state_covariance.size2(), 3u)
      << "current_state_covariance columns were not sized by the dimensions.";
}

/**
 * @test A TransitionState built from a dimensions object must size both
 * matrices from that object.
 */
TEST(FilterStatesValidationTest, transitionStateSizesBuffersFromDimensions) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  const TransitionState transition_state(dimensions);
  EXPECT_EQ(transition_state.transition_matrix.size1(), 3u)
      << "transition_matrix rows were not sized by the dimensions.";
  EXPECT_EQ(transition_state.transition_covariance.size2(), 3u)
      << "transition_covariance columns were not sized by the dimensions.";
}

/**
 * @test The std::vector overload of setTransitionCovariance must reject a
 * source with more rows than the target.
 */
TEST(
    FilterStatesValidationTest, setTransitionCovarianceVectorRejectsRowMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const std::vector<std::vector<double>> oversized(
      4, std::vector<double>(3, 1.0)
  );
  EXPECT_THROW(
      kca_states.setTransitionCovariance(oversized), filter_shape_mismatch
  ) << "setTransitionCovariance accepted a 4x3 source for a 3x3 target.";
}

/**
 * @test The std::vector overload of setCurrentStateCovariance must reject a
 * source with more columns than the target.
 */
TEST(
    FilterStatesValidationTest,
    setCurrentStateCovarianceVectorRejectsColumnMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const std::vector<std::vector<double>> oversized(
      3, std::vector<double>(4, 1.0)
  );
  EXPECT_THROW(
      kca_states.setCurrentStateCovariance(oversized), filter_shape_mismatch
  ) << "setCurrentStateCovariance accepted a 3x4 source for a 3x3 target.";
}

/**
 * @test The std::vector overload of setObservationMatrix must reject a source
 * with more rows than the 1x3 target.
 */
TEST(FilterStatesValidationTest, setObservationMatrixVectorRejectsRowMismatch) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const std::vector<std::vector<double>> oversized(
      2, std::vector<double>(3, 1.0)
  );
  EXPECT_THROW(
      kca_states.setObservationMatrix(oversized), filter_shape_mismatch
  ) << "setObservationMatrix accepted a 2x3 source for a 1x3 target.";
}

/**
 * @test The std::vector overload of setPredictedObservationCovariance must
 * reject a source larger than the 1x1 target.
 */
TEST(
    FilterStatesValidationTest,
    setPredictedObservationCovarianceVectorRejectsShapeMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const std::vector<std::vector<double>> oversized(
      2, std::vector<double>(2, 1.0)
  );
  EXPECT_THROW(
      kca_states.setPredictedObservationCovariance(oversized),
      filter_shape_mismatch
  ) << "setPredictedObservationCovariance accepted a 2x2 source for a 1x1 "
       "target.";
}

/**
 * @test The std::vector overload of setPredictedObservationMean must reject a
 * source longer than the length-1 target.
 */
TEST(
    FilterStatesValidationTest,
    setPredictedObservationMeanVectorRejectsLengthMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const std::vector<double> oversized(2, 1.0);
  EXPECT_THROW(
      kca_states.setPredictedObservationMean(oversized), filter_shape_mismatch
  ) << "setPredictedObservationMean accepted a length-2 source for a length-1 "
       "target.";
}

/**
 * @test The std::vector overload of setPredictedStateCovariance must reject a
 * source with more rows than the target.
 */
TEST(
    FilterStatesValidationTest,
    setPredictedStateCovarianceVectorRejectsRowMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const std::vector<std::vector<double>> oversized(
      4, std::vector<double>(3, 1.0)
  );
  EXPECT_THROW(
      kca_states.setPredictedStateCovariance(oversized), filter_shape_mismatch
  ) << "setPredictedStateCovariance accepted a 4x3 source for a 3x3 target.";
}

/**
 * @test The std::vector overload of setPredictedStateMean must reject a source
 * longer than the target.
 */
TEST(
    FilterStatesValidationTest, setPredictedStateMeanVectorRejectsLengthMismatch
) {
  const FilterSystemDimensions dimensions(3, 3, 3, 1, 3, 1, 1, 0.0);
  KcaStates kca_states(dimensions);
  const std::vector<double> oversized(4, 1.0);
  EXPECT_THROW(
      kca_states.setPredictedStateMean(oversized), filter_shape_mismatch
  ) << "setPredictedStateMean accepted a length-4 source for a length-3 "
       "target.";
}
