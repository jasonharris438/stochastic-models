#include "stochastic_models/kalman_filter/states.h"
#include "stochastic_models/kalman_filter/type_conversion.h"
#include "stochastic_models/numeric_utils/linalg.h"

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <cmath>
#include <gtest/gtest.h>
#include <vector>

using namespace boost::numeric::ublas;
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

  const std::vector<std::vector<double>> expected_result{
      {{0.013744, 0.0, 0.0}, {0.0, 0.001, 0.0}, {0.0, 0.0, 0.001}}
  };
  EXPECT_EQ(result_vector, expected_result) << "The covariance matrix is not "
                                               "calculated correctly.";
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
  EXPECT_EQ(result.size(), boost_matrix.size2())
      << "The mean vector is of incorrect dimension.";

  const std::vector<double> result_vector{result(0), result(1), result(2)};
  const std::vector<double> expected_result{1.3309855249349998, 0, 0};
  EXPECT_EQ(result_vector, expected_result)
      << "The mean vector is not calculated correctly.";
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
  const std::vector<double> expected_result{1.330986};
  EXPECT_EQ(result.size(), expected_result.size())
      << "The mean vector is of incorrect dimension.";

  const std::vector<double> result_vector{result(0)};
  EXPECT_EQ(result_vector, expected_result)
      << "The mean vector is not calculated correctly.";
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

  // Convert result to a std::vector from boost types.
  const std::vector<std::vector<double>> result_vec =
      copy_matrix_elements_to_vector(result);

  // Expected result from above operations.
  const std::vector<std::vector<double>> expected_result{{0.0137912692000676}};

  // Assert shape of result.
  EXPECT_EQ(result.size1(), expected_result.size())
      << "The covariance matrix is of incorrect row dimension.";
  EXPECT_EQ(result.size2(), expected_result.at(0).size())
      << "The covariance matrix is of incorrect column dimension.";

  // Assert value(s) of result.
  const std::vector<double> result_vector{result(0)};
  EXPECT_EQ(result_vec, expected_result)
      << "The covariance vector is not calculated correctly.";
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

  // Generate predicted observation covariance from variables created.
  const BoostMartixInverter matrix_inverter;
  const matrix<double> result = predicted_observation.calculateKalmanGain(
      predicted_state_covariance, predicted_observation_covariance,
      matrix_inverter
  );

  // Convert result to a std::vector from boost types.
  const std::vector<std::vector<double>> result_vec =
      copy_matrix_elements_to_vector(result);

  // Expected result from above operations.
  const std::vector<std::vector<double>> expected_result{
      {0.99654139808725539}, {0}, {0}
  };

  // Assert shape of result.
  EXPECT_EQ(result.size1(), expected_result.size())
      << "The covariance matrix is of incorrect row dimension.";
  EXPECT_EQ(result.size2(), expected_result.at(0).size())
      << "The covariance matrix is of incorrect column dimension.";

  // Assert value(s) of result.
  const std::vector<double> result_vector{result(0)};
  EXPECT_EQ(result_vec, expected_result)
      << "The covariance vector is not calculated correctly.";
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
  const std::vector<double> expected_result{1.3108695137597879, 0, 0};

  EXPECT_EQ(result.size(), expected_result.size())
      << "The mean vector is of incorrect dimension.";

  const std::vector<double> result_vector{result(0), result(1), result(2)};
  EXPECT_EQ(result_vector, expected_result)
      << "The mean vector is not calculated correctly.";
}

/**
 * @brief Test that the CurrentState.calculateMean method returns the
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
  const std::vector<std::vector<double>> expected_result{
      {4.7105773280001206e-05, 0, 0}, {0, 0.001, 0}, {0, 0, 0.001}
  };

  EXPECT_EQ(result_vec.size(), expected_result.size())
      << "The covariance matrix is of incorrect dimension.";

  EXPECT_EQ(result_vec, expected_result)
      << "The covariance matrix is not calculated correctly.";
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
  std::vector<double> predicted_state_mean_vector(3);
  std::copy(
      predicted_state_mean.begin(), predicted_state_mean.end(),
      predicted_state_mean_vector.begin()
  );
  const std::vector<double> expected_predicted_state_mean{
      10.301048359829961, 0.0, 0.0
  };
  EXPECT_EQ(predicted_state_mean_vector, expected_predicted_state_mean)
      << "The predicted state mean vector was set with invalid or "
         "inconsistent values.";

  // Predicted state covariance.
  const std::vector<std::vector<double>> predicted_state_covariance =
      copy_matrix_elements_to_vector(kca_states.getPredictedStateCovariance());
  const std::vector<std::vector<double>> expected_predicted_state_covariance{
      {0.12695229227341848, 0, 0}, {0, 0.001, 0}, {0, 0, 0.001}
  };
  EXPECT_EQ(predicted_state_covariance, expected_predicted_state_covariance)
      << "The predicted state covariance was set with invalid or "
         "inconsistent values.";
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
  std::vector<double> predicted_observation_mean{10.301, 0, 0};

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
  std::vector<double> predicted_observation_mean_vector(3);
  std::copy(
      kca_states.getPredictedObservationMean().begin(),
      kca_states.getPredictedObservationMean().end(),
      predicted_observation_mean_vector.begin()
  );
  const std::vector<double> expected_predicted_observation_mean{
      10.301048359829961, 0, 0
  };
  EXPECT_EQ(
      predicted_observation_mean_vector, expected_predicted_observation_mean
  ) << "The predicted observation mean vector was set with invalid or "
       "inconsistent values.";

  // Predicted observation covariance.
  const std::vector<std::vector<double>>
      predicted_observation_covariance_vector = copy_matrix_elements_to_vector(
          kca_states.getPredictedObservationCovariance()
      );
  const std::vector<std::vector<double>>
      expected_predicted_observation_covariance{{0.13695229227341849}};
  EXPECT_EQ(
      predicted_observation_covariance_vector,
      expected_predicted_observation_covariance
  ) << "The predicted observation covariance was set with invalid or "
       "inconsistent values.";

  // Current state mean.
  const vector<double> new_current_state_mean =
      kca_states.getCurrentStateMean();
  std::vector<double> new_current_state_mean_vector(3);
  std::copy(
      new_current_state_mean.begin(), new_current_state_mean.end(),
      new_current_state_mean_vector.begin()
  );
  const std::vector<double> expected_current_state_mean{10.3000765492722, 0, 0};
  EXPECT_EQ(new_current_state_mean_vector, expected_current_state_mean)
      << "The current state mean vector was set with invalid or "
         "inconsistent values.";

  // Current state covariance.
  const std::vector<std::vector<double>> new_current_state_covariance =
      copy_matrix_elements_to_vector(kca_states.getCurrentStateCovariance());
  const std::vector<std::vector<double>> expected_current_state_covariance{
      {0.0092698187205194488, 0, 0}, {0, 0.001, 0}, {0, 0, 0.001}
  };
  EXPECT_EQ(new_current_state_covariance, expected_current_state_covariance)
      << "The current state covariance was set with invalid or "
         "inconsistent values.";
}
