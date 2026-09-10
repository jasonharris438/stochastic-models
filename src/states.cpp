#include "stochastic_models/kalman_filter/states.h"

#include "stochastic_models/kalman_filter/states_exceptions.h"
#include "stochastic_models/kalman_filter/type_conversion.h"
#include "stochastic_models/numeric_utils/helpers.h"
#include "stochastic_models/numeric_utils/linalg.h"

#include <algorithm>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <cmath>
#include <cstddef>
#include <string>
#include <string_view>

// Just for this module as we do not introduce any other namespaces.
using namespace boost::numeric::ublas;

namespace {

  void requireDimensionInRange(std::string_view name, std::size_t value) {
    if (value == 0 || value > FilterSystemDimensions::max_dimension) {
      throw invalid_filter_dimensions(
          "Dimension '" + std::string{name} + "' must be in [1, " +
          std::to_string(FilterSystemDimensions::max_dimension) + "]; got " +
          std::to_string(value) + "."
      );
    }
  }

  void requireDimensionsEqual(
      std::string_view name,
      std::size_t value,
      std::string_view reference_name,
      std::size_t reference
  ) {
    if (value != reference) {
      throw invalid_filter_dimensions(
          "Dimension '" + std::string{name} + "' (" + std::to_string(value) +
          ") must equal '" + std::string{reference_name} + "' (" +
          std::to_string(reference) + ")."
      );
    }
  }

  void requireScalarObservation(const FilterSystemDimensions& dimensions) {
    requireDimensionsEqual(
        "observation_matrix_rows", dimensions.getObservationMatrixRows(),
        "observation_dimension", KcaStates::observation_dimension
    );
  }

  void requireLength(
      std::string_view name,
      std::size_t source_length,
      std::size_t target_length
  ) {
    if (source_length != target_length) {
      throw filter_shape_mismatch(
          "Source length " + std::to_string(source_length) + " for '" +
          std::string{name} + "' does not match target length " +
          std::to_string(target_length) + "."
      );
    }
  }

  void assignChecked(
      std::string_view name,
      const vector<double>& source,
      vector<double>& target
  ) {
    requireLength(name, source.size(), target.size());
    target.assign(source);
  }

  void assignChecked(
      std::string_view name,
      const std::vector<double>& source,
      vector<double>& target
  ) {
    requireLength(name, source.size(), target.size());
    std::copy(source.begin(), source.end(), target.begin());
  }

  void requireShape(
      std::string_view name,
      std::size_t source_rows,
      std::size_t source_columns,
      const matrix<double>& target
  ) {
    if (source_rows != target.size1() || source_columns != target.size2()) {
      throw filter_shape_mismatch(
          "Source shape " + std::to_string(source_rows) + "x" +
          std::to_string(source_columns) + " for '" + std::string{name} +
          "' does not match target shape " + std::to_string(target.size1()) +
          "x" + std::to_string(target.size2()) + "."
      );
    }
  }

  void assignChecked(
      std::string_view name,
      const matrix<double>& source,
      matrix<double>& target
  ) {
    requireShape(name, source.size1(), source.size2(), target);
    target.assign(source);
  }

  void assignChecked(
      std::string_view name,
      const std::vector<std::vector<double>>& source,
      matrix<double>& target
  ) {
    const std::size_t source_columns =
        source.empty() ? 0 : source.front().size();
    requireShape(name, source.size(), source_columns, target);
    for (const std::vector<double>& source_row : source) {
      requireLength(name, source_row.size(), target.size2());
    }
    for (std::size_t i{0}; i < target.size1(); i++) {
      std::copy(source[i].begin(), source[i].end(), row(target, i).begin());
    }
  }

} // namespace

// Prior predicted state class functionality implementation.
PredictedState::PredictedState(matrix<double> transition_matrix)
    : transition_matrix(transition_matrix) {}

void PredictedState::setTransitionMatrix(
    const matrix<double>& transition_matrix
) {
  this->transition_matrix = transition_matrix;
}
const vector<double>
PredictedState::calculateMean(const vector<double>& current_state_mean) const {
  const vector<double> new_mean = prod(transition_matrix, current_state_mean);
  return new_mean;
}
const matrix<double> PredictedState::calculateCovariance(
    const matrix<double>& current_state_covariance,
    const matrix<double>& transition_covariance
) const {
  const matrix<double> current_state_transition_matrix =
      prod(current_state_covariance, trans(transition_matrix));
  const matrix<double> new_covariance =
      prod(transition_matrix, current_state_transition_matrix) +
      transition_covariance;
  return new_covariance;
}

// Prior predicted observation class functionality implementation.
PredictedObservation::PredictedObservation(
    matrix<double> observation_matrix, double observation_offset
)
    : observation_matrix(observation_matrix),
      observation_offset(observation_offset) {}
const matrix<double>& PredictedObservation::getObservationMatrix() const {
  return observation_matrix;
}
void PredictedObservation::setObservationMatrix(
    const matrix<double>& observation_matrix
) {
  this->observation_matrix = observation_matrix;
}
void PredictedObservation::setObservationOffset(
    const double& observation_offset
) {
  this->observation_offset = observation_offset;
}
const vector<double> PredictedObservation::calculateMean(
    const vector<double>& predicted_state_mean
) const {
  vector<double> new_mean = prod(observation_matrix, predicted_state_mean);
  add_scalar_to_vector(new_mean, observation_offset);
  return new_mean;
}
const matrix<double> PredictedObservation::calculateCovariance(
    const matrix<double>& predicted_state_covariance,
    const double& innovation_sigma
) const {
  // Result matrix components.
  matrix<double> observation_transposed = trans(observation_matrix);
  matrix<double> inner_product =
      prod(predicted_state_covariance, observation_transposed);
  matrix<double> result = prod(observation_matrix, inner_product);

  // Add the squared innovation sigma to the result.
  const double innovation_sigma_squared = std::pow(innovation_sigma, 2);
  add_scalar_to_matrix(result, innovation_sigma_squared);

  return result;
}

// Posterior predicted current state class functionality implementation.
const vector<double> CurrentState::calculateMean(
    const vector<double>& predicted_state_mean,
    const vector<double>& kalman_gain,
    const double& innovation
) const {
  const vector<double> inner_product = kalman_gain * innovation;
  vector<double> result = predicted_state_mean + inner_product;
  return result;
}
const matrix<double> CurrentState::calculateCovariance(
    const matrix<double>& predicted_state_covariance,
    const matrix<double>& observation_matrix,
    const matrix<double>& kalman_gain
) const {
  const matrix<double> inner_matrix_product =
      prod(observation_matrix, predicted_state_covariance);
  const matrix<double> filter_covariance_product =
      prod(kalman_gain, inner_matrix_product);
  const matrix<double> result =
      predicted_state_covariance - filter_covariance_product;
  return result;
}

// Prior state data class / struct implementation.
PriorState::PriorState(
    std::size_t state_mean_dimension,
    std::size_t state_covariance_rows,
    std::size_t state_covariance_columns,
    std::size_t observation_matrix_rows,
    std::size_t observation_matrix_columns,
    std::size_t observation_covariance_rows,
    std::size_t observation_covariance_columns,
    const double& observation_offset
)
    : predicted_observation_mean(vector<double>(observation_matrix_rows)),
      predicted_state_mean(vector<double>(state_mean_dimension)),
      predicted_observation_covariance(
          matrix<double>(
              observation_covariance_rows, observation_covariance_columns
          )
      ),
      predicted_state_covariance(
          matrix<double>(state_covariance_rows, state_covariance_columns)
      ),
      observation_matrix(
          matrix<double>(observation_matrix_rows, observation_matrix_columns)
      ),
      observation_offset(observation_offset) {}

// Posterior state data class / struct implementation.
PosteriorState::PosteriorState(
    std::size_t state_mean_dimension,
    std::size_t state_covariance_rows,
    std::size_t state_covariance_columns
)
    : current_state_mean(vector<double>(state_mean_dimension)),
      current_state_covariance(
          matrix<double>(state_covariance_rows, state_covariance_columns)
      ) {}

// Transition state data class / struct implementation
TransitionState::TransitionState(
    std::size_t state_covariance_rows, std::size_t state_covariance_columns
)
    : transition_matrix(
          matrix<double>(state_covariance_rows, state_covariance_columns)
      ),
      transition_covariance(
          matrix<double>(state_covariance_rows, state_covariance_columns)
      ) {}

// Filter boolean state data class / struct implementation
FilterState::FilterState() : initialised(false), priors_set(false) {}

// Type that contains the validated dimensions of a Kalman Filter system.
FilterSystemDimensions::FilterSystemDimensions(
    std::size_t state_mean_dimension,
    std::size_t state_covariance_rows,
    std::size_t state_covariance_columns,
    std::size_t observation_matrix_rows,
    std::size_t observation_matrix_columns,
    std::size_t observation_covariance_rows,
    std::size_t observation_covariance_columns,
    double observation_offset
)
    : state_mean_dimension(state_mean_dimension),
      state_covariance_rows(state_covariance_rows),
      state_covariance_columns(state_covariance_columns),
      observation_matrix_rows(observation_matrix_rows),
      observation_matrix_columns(observation_matrix_columns),
      observation_covariance_rows(observation_covariance_rows),
      observation_covariance_columns(observation_covariance_columns),
      observation_offset(observation_offset) {
  requireDimensionInRange("state_mean_dimension", state_mean_dimension);
  requireDimensionInRange("state_covariance_rows", state_covariance_rows);
  requireDimensionInRange("state_covariance_columns", state_covariance_columns);
  requireDimensionInRange("observation_matrix_rows", observation_matrix_rows);
  requireDimensionInRange(
      "observation_matrix_columns", observation_matrix_columns
  );
  requireDimensionInRange(
      "observation_covariance_rows", observation_covariance_rows
  );
  requireDimensionInRange(
      "observation_covariance_columns", observation_covariance_columns
  );
  requireDimensionsEqual(
      "state_covariance_rows", state_covariance_rows, "state_mean_dimension",
      state_mean_dimension
  );
  requireDimensionsEqual(
      "state_covariance_columns", state_covariance_columns,
      "state_mean_dimension", state_mean_dimension
  );
  requireDimensionsEqual(
      "observation_matrix_columns", observation_matrix_columns,
      "state_mean_dimension", state_mean_dimension
  );
  requireDimensionsEqual(
      "observation_covariance_rows", observation_covariance_rows,
      "observation_matrix_rows", observation_matrix_rows
  );
  requireDimensionsEqual(
      "observation_covariance_columns", observation_covariance_columns,
      "observation_matrix_rows", observation_matrix_rows
  );
}
std::size_t FilterSystemDimensions::getStateMeanDimension() const noexcept {
  return state_mean_dimension;
}
std::size_t FilterSystemDimensions::getStateCovarianceRows() const noexcept {
  return state_covariance_rows;
}
std::size_t FilterSystemDimensions::getStateCovarianceColumns() const noexcept {
  return state_covariance_columns;
}
std::size_t FilterSystemDimensions::getObservationMatrixRows() const noexcept {
  return observation_matrix_rows;
}
std::size_t
FilterSystemDimensions::getObservationMatrixColumns() const noexcept {
  return observation_matrix_columns;
}
std::size_t
FilterSystemDimensions::getObservationCovarianceRows() const noexcept {
  return observation_covariance_rows;
}
std::size_t
FilterSystemDimensions::getObservationCovarianceColumns() const noexcept {
  return observation_covariance_columns;
}
double FilterSystemDimensions::getObservationOffset() const noexcept {
  return observation_offset;
}

// State handler for the KCA implementation.
KcaStates::KcaStates(const FilterSystemDimensions& dimensions)
    : prior_state(
          dimensions.getStateMeanDimension(),
          dimensions.getStateCovarianceRows(),
          dimensions.getStateCovarianceColumns(),
          dimensions.getObservationMatrixRows(),
          dimensions.getObservationMatrixColumns(),
          dimensions.getObservationCovarianceRows(),
          dimensions.getObservationCovarianceColumns(),
          dimensions.getObservationOffset()
      ),
      posterior_state(
          dimensions.getStateMeanDimension(),
          dimensions.getStateCovarianceRows(),
          dimensions.getStateCovarianceColumns()
      ),
      transition_state(
          dimensions.getStateCovarianceRows(),
          dimensions.getStateCovarianceColumns()
      ) {
  requireScalarObservation(dimensions);
}

void KcaStates::setInitialState(
    const std::vector<double>& data_series, const double& h, const double& q
) {
  check_minimum_observations(data_series, 1, "KCA initialisation");

  // Initial transition state.
  std::vector<std::vector<double>> transition_matrix_as_vectors{
      {1.0, h, 0.5 * std::pow(h, 2)}, {0.0, 1.0, h}, {0.0, 0.0, 1.0}
  };
  std::vector<std::vector<double>> transition_covariance_as_vectors{
      {q * std::pow(h, 5) / 20.0, q * std::pow(h, 4) / 8.0,
       q * std::pow(h, 3) / 6.0},
      {q * std::pow(h, 4) / 8.0, q * std::pow(h, 3) / 3.0,
       q * std::pow(h, 2) / 2.0},
      {q * std::pow(h, 3) / 6.0, q * std::pow(h, 2) / 2.0, q * h}
  };

  // Initial current state.
  std::vector<double> current_state_mean_as_vector{
      data_series.back(), 0.0, 0.0
  };
  std::vector<std::vector<double>> current_state_covariance_as_vectors{
      {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}
  };

  // Initial observation state.
  std::vector<std::vector<double>> observation_matrix_as_vectors{
      {1.0, 0.0, 0.0}
  };
  const double observation_offset = 0.0;

  // Move the current state and transition state to the target matrices.
  setCurrentStateMean(current_state_mean_as_vector);
  setCurrentStateCovariance(current_state_covariance_as_vectors);

  // Move the the transition and transition covariance matrices to the
  // target objects.
  setTransitionMatrix(transition_matrix_as_vectors);
  setTransitionCovariance(transition_covariance_as_vectors);

  // Move the the observation matrix to the target matrix, and copy the
  // observation offset.
  setObservationMatrix(observation_matrix_as_vectors);
  setObservationOffset(observation_offset);

  // We are now fully initialised.
  setInitialized();
}
void KcaStates::updatePredictedState() {
  if (!isInitialised()) {
    throw filter_uninitialised(
        "The KCA kalman filter has not been initialised."
    );
  }

  const PredictedState predicted_state{getTransitionMatrix()};

  const vector<double> predicted_state_mean =
      predicted_state.calculateMean(getCurrentStateMean());
  const matrix<double> predicted_state_covariance =
      predicted_state.calculateCovariance(
          getCurrentStateCovariance(), getTransitionCovariance()
      );

  setPredictedStateMean(predicted_state_mean);
  setPredictedStateCovariance(predicted_state_covariance);

  // Set the priors to true after the predicted state has been updated.
  setPriorsTrue();
}
void KcaStates::updateCurrentState(
    const double& observation, const double& innovation_sigma
) {
  if (!isInitialised()) {
    throw filter_uninitialised(
        "The KCA kalman filter has not been initialised."
    );
  }
  if (!arePriorsValid()) {
    throw filter_invalid_operation(
        "The KCA kalman filter priors must be set to valid state "
        "before calling updateCurrentState."
    );
  }

  const PredictedObservation predicted_observation{
      getObservationMatrix(), getObservationOffset()
  };

  const vector<double> predicted_observation_mean =
      predicted_observation.calculateMean(getPredictedStateMean());
  const matrix<double> predicted_observation_covariance =
      predicted_observation.calculateCovariance(
          getPredictedStateCovariance(), innovation_sigma
      );

  const BoostMatrixInverter matrix_inverter;
  const matrix<double> kalman_gain = predicted_observation.calculateKalmanGain(
      getPredictedStateCovariance(), predicted_observation_covariance,
      matrix_inverter
  );
  const double innovation = observation - predicted_observation_mean(0);

  vector<double> kalman_gain_vector = vector<double>(kalman_gain.size1());
  std::copy(
      column(kalman_gain, 0).begin(), column(kalman_gain, 0).end(),
      kalman_gain_vector.begin()
  );

  const CurrentState current_state;
  const vector<double> current_state_mean = current_state.calculateMean(
      getPredictedStateMean(), kalman_gain_vector, innovation
  );
  const matrix<double> current_state_covariance =
      current_state.calculateCovariance(
          getPredictedStateCovariance(),
          predicted_observation.getObservationMatrix(), kalman_gain
      );

  // Calculate everything first to ensure that state is not half-set when
  // an exception is thrown.
  setPredictedObservationMean(predicted_observation_mean);
  setPredictedObservationCovariance(predicted_observation_covariance);
  setCurrentStateMean(current_state_mean);
  setCurrentStateCovariance(current_state_covariance);

  // Priors are now in an invalid state for a further posterior update.
  setPriorsFalse();
}
const vector<double>& KcaStates::getCurrentStateMean() const {
  return posterior_state.current_state_mean;
}
const std::vector<double> KcaStates::getCurrentStateMeanVector() const {
  std::vector<double> current_state_mean_vector =
      std::vector<double>(posterior_state.current_state_mean.size());
  std::copy(
      posterior_state.current_state_mean.begin(),
      posterior_state.current_state_mean.end(),
      current_state_mean_vector.begin()
  );

  return current_state_mean_vector;
}
const matrix<double>& KcaStates::getCurrentStateCovariance() const {
  return posterior_state.current_state_covariance;
}
const matrix<double>& KcaStates::getObservationMatrix() const {
  return prior_state.observation_matrix;
}
const double& KcaStates::getObservationOffset() const {
  return prior_state.observation_offset;
}
const matrix<double>& KcaStates::getPredictedObservationCovariance() const {
  return prior_state.predicted_observation_covariance;
}
const vector<double>& KcaStates::getPredictedObservationMean() const {
  return prior_state.predicted_observation_mean;
}
const matrix<double>& KcaStates::getPredictedStateCovariance() const {
  return prior_state.predicted_state_covariance;
}
const vector<double>& KcaStates::getPredictedStateMean() const {
  return prior_state.predicted_state_mean;
}
const matrix<double>& KcaStates::getTransitionCovariance() const {
  return transition_state.transition_covariance;
}
const matrix<double>& KcaStates::getTransitionMatrix() const {
  return transition_state.transition_matrix;
}
const bool& KcaStates::isInitialised() const {
  return filter_state.initialised;
}
const bool& KcaStates::arePriorsValid() const {
  return filter_state.priors_set;
}
void KcaStates::setInitialized() {
  filter_state.initialised = true;
}
void KcaStates::setPriorsTrue() {
  filter_state.priors_set = true;
}
void KcaStates::setPriorsFalse() {
  filter_state.priors_set = false;
}
void KcaStates::setCurrentStateMean(const vector<double>& current_state_mean) {
  assignChecked(
      "current_state_mean", current_state_mean,
      posterior_state.current_state_mean
  );
}
void KcaStates::setCurrentStateMean(
    const std::vector<double>& current_state_mean
) {
  assignChecked(
      "current_state_mean", current_state_mean,
      posterior_state.current_state_mean
  );
}
void KcaStates::setCurrentStateCovariance(
    const matrix<double>& current_state_covariance
) {
  assignChecked(
      "current_state_covariance", current_state_covariance,
      posterior_state.current_state_covariance
  );
}
void KcaStates::setCurrentStateCovariance(
    const std::vector<std::vector<double>>& current_state_covariance
) {
  assignChecked(
      "current_state_covariance", current_state_covariance,
      posterior_state.current_state_covariance
  );
}
void KcaStates::setObservationMatrix(const matrix<double>& observation_matrix) {
  assignChecked(
      "observation_matrix", observation_matrix, prior_state.observation_matrix
  );
}
void KcaStates::setObservationMatrix(
    const std::vector<std::vector<double>>& observation_matrix
) {
  assignChecked(
      "observation_matrix", observation_matrix, prior_state.observation_matrix
  );
}
void KcaStates::setObservationOffset(const double& observation_offset) {
  prior_state.observation_offset = observation_offset;
}
void KcaStates::setPredictedObservationCovariance(
    const matrix<double>& predicted_observation_covariance
) {
  assignChecked(
      "predicted_observation_covariance", predicted_observation_covariance,
      prior_state.predicted_observation_covariance
  );
}
void KcaStates::setPredictedObservationCovariance(
    const std::vector<std::vector<double>>& predicted_observation_covariance
) {
  assignChecked(
      "predicted_observation_covariance", predicted_observation_covariance,
      prior_state.predicted_observation_covariance
  );
}
void KcaStates::setPredictedObservationMean(
    const vector<double>& predicted_observation_mean
) {
  assignChecked(
      "predicted_observation_mean", predicted_observation_mean,
      prior_state.predicted_observation_mean
  );
}
void KcaStates::setPredictedObservationMean(
    const std::vector<double>& predicted_observation_mean
) {
  assignChecked(
      "predicted_observation_mean", predicted_observation_mean,
      prior_state.predicted_observation_mean
  );
}
void KcaStates::setPredictedStateCovariance(
    const matrix<double>& predicted_state_covariance
) {
  assignChecked(
      "predicted_state_covariance", predicted_state_covariance,
      prior_state.predicted_state_covariance
  );
}
void KcaStates::setPredictedStateCovariance(
    const std::vector<std::vector<double>>& predicted_state_covariance
) {
  assignChecked(
      "predicted_state_covariance", predicted_state_covariance,
      prior_state.predicted_state_covariance
  );
}
void KcaStates::setPredictedStateMean(
    const vector<double>& predicted_state_mean
) {
  assignChecked(
      "predicted_state_mean", predicted_state_mean,
      prior_state.predicted_state_mean
  );
}
void KcaStates::setPredictedStateMean(
    const std::vector<double>& predicted_state_mean
) {
  assignChecked(
      "predicted_state_mean", predicted_state_mean,
      prior_state.predicted_state_mean
  );
}
void KcaStates::setTransitionCovariance(
    const matrix<double>& transition_covariance
) {
  assignChecked(
      "transition_covariance", transition_covariance,
      transition_state.transition_covariance
  );
}
void KcaStates::setTransitionCovariance(
    const std::vector<std::vector<double>>& transition_covariance
) {
  assignChecked(
      "transition_covariance", transition_covariance,
      transition_state.transition_covariance
  );
}
void KcaStates::setTransitionMatrix(const matrix<double>& transition_matrix) {
  assignChecked(
      "transition_matrix", transition_matrix, transition_state.transition_matrix
  );
}
void KcaStates::setTransitionMatrix(
    const std::vector<std::vector<double>>& transition_matrix
) {
  assignChecked(
      "transition_matrix", transition_matrix, transition_state.transition_matrix
  );
}
