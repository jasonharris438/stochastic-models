#include "stochastic_models/kalman_filter/states.h"

#include "stochastic_models/kalman_filter/states_exceptions.h"
#include "stochastic_models/kalman_filter/type_conversion.h"
#include "stochastic_models/numeric_utils/linalg.h"

#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <cmath>
#include <utility>

// Just for this module as we do not introduce any other namespaces.
using namespace boost::numeric::ublas;

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
    const int& state_mean_dimension,
    const int& state_covariance_rows,
    const int& state_covariance_columns,
    const int& observation_matrix_rows,
    const int& observation_matrix_columns,
    const int& observation_covariance_rows,
    const int& observation_covariance_columns,
    const double& observation_offset
)
    : predicted_observation_mean(vector<double>(state_mean_dimension)),
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
    const int& state_mean_dimension,
    const int& state_covariance_rows,
    const int& state_covariance_columns
)
    : current_state_mean(vector<double>(state_mean_dimension)),
      current_state_covariance(
          matrix<double>(state_covariance_rows, state_covariance_columns)
      ) {}

// Transition state data class / struct implementation
TransitionState::TransitionState(
    const int& state_covariance_rows, const int& state_covariance_columns
)
    : transition_matrix(
          matrix<double>(state_covariance_rows, state_covariance_columns)
      ),
      transition_covariance(
          matrix<double>(state_covariance_rows, state_covariance_columns)
      ) {}

// Filter boolean state data class / struct implementation
FilterState::FilterState() : initialised(false), priors_set(false) {}

// General SDE state handler for managing stochastic process that governs series
// being analysed by the KCA.
FilterGeneralSde::FilterGeneralSde(
    const double& mu,
    const double& sigma,
    const double& conditional_variance,
    const double& mu_numerator,
    const double& mu_denominator,
    const double& sigma_kernel_value,
    const u_int32_t& n_observations
)
    : mu(mu), sigma(sigma), conditional_variance(conditional_variance),
      mu_numerator(mu_numerator), mu_denominator(mu_denominator),
      sigma_kernel_value(sigma_kernel_value), n_observations(n_observations) {}
FilterGeneralSde::FilterGeneralSde()
    : FilterGeneralSde(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0) {}
const double& FilterGeneralSde::getMu() const {
  return mu;
}
const double& FilterGeneralSde::getSigma() const {
  return sigma;
}
const double& FilterGeneralSde::getConditionalVariance() const {
  return conditional_variance;
}
void FilterGeneralSde::initializeLikelihoodState(
    const std::vector<double>& data_series
) {
  const std::unordered_map<std::string, const double> likelihood_values =
      likelihood.calculate(data_series);
  mu = likelihood_values.at("mu");
  sigma = likelihood_values.at("sigma");
  conditional_variance = likelihood.calculateConditionalVariance(sigma, mu);

  u_int32_t n_obs = static_cast<u_int32_t>(data_series.size());

  n_observations = n_obs;
  mu_numerator = likelihood.calculateLeadLagInnerProduct(data_series);
  mu_denominator = likelihood.calculateLagSquared(data_series);
  sigma_kernel_value = likelihood.calculateSigmaKernel(n_obs, sigma);
}

// Type that contains internal dimensions of a Kalman Filter system.
FilterSystemDimensions::FilterSystemDimensions()
    : state_mean_dimension(0), state_covariance_rows(0),
      state_covariance_columns(0), observation_matrix_rows(0),
      observation_matrix_columns(0), observation_covariance_rows(0),
      observation_covariance_columns(0), observation_offset(0.0) {}
FilterSystemDimensions::FilterSystemDimensions(
    int state_mean_dimension,
    int state_covariance_rows,
    int state_covariance_columns,
    int observation_matrix_rows,
    int observation_matrix_columns,
    int observation_covariance_rows,
    int observation_covariance_columns,
    double observation_offset
)
    : state_mean_dimension(state_mean_dimension),
      state_covariance_rows(state_covariance_rows),
      state_covariance_columns(state_covariance_columns),
      observation_matrix_rows(observation_matrix_rows),
      observation_matrix_columns(observation_matrix_columns),
      observation_covariance_rows(observation_covariance_rows),
      observation_covariance_columns(observation_covariance_columns),
      observation_offset(observation_offset) {}

// State handler for the KCA implementation.
KcaStates::KcaStates(const FilterSystemDimensions& dimensions)
    : prior_state(
          dimensions.state_mean_dimension,
          dimensions.state_covariance_rows,
          dimensions.state_covariance_columns,
          dimensions.observation_matrix_rows,
          dimensions.observation_matrix_columns,
          dimensions.observation_covariance_rows,
          dimensions.observation_covariance_columns,
          dimensions.observation_offset
      ),
      posterior_state(
          dimensions.state_mean_dimension,
          dimensions.state_covariance_rows,
          dimensions.state_covariance_columns
      ),
      transition_state(
          dimensions.state_covariance_rows, dimensions.state_covariance_columns
      ) {}

void KcaStates::move_std_vectors_to_matrix(
    std::vector<std::vector<double>>&& matrix_as_vectors, matrix<double>& target
) {
  u_int32_t rows = static_cast<u_int32_t>(matrix_as_vectors.size());

  // Move the vectors into the target matrix.
  for (u_int32_t i{0}; i < rows; i++) {
    std::move(
        matrix_as_vectors.at(i).begin(), matrix_as_vectors.at(i).end(),
        row(target, i).begin()
    );
  }
}
void KcaStates::move_std_vector_to_vector(
    std::vector<double>&& vector_as_vector, vector<double>& target
) {
  // Move the vector into the target vector.
  std::move(vector_as_vector.begin(), vector_as_vector.end(), target.begin());
}
void KcaStates::setInitialState(
    const std::vector<double>& data_series, const double& h, const double& q
) {
  filter_sde.initializeLikelihoodState(data_series);

  // Initial transition state.
  std::vector<std::vector<double>> transition_matrix_as_vectors{
      {std::exp(filter_sde.getMu()), h, 0.5 * std::pow(h, 2)},
      {0.0, 1.0, h},
      {0.0, 0.0, 1.0}
  };
  std::vector<std::vector<double>> transition_covariance_as_vectors{
      {filter_sde.getConditionalVariance(), 0.0, 0.0},
      {0.0, q, 0.0},
      {0.0, 0.0, q}
  };

  // Initial current state.
  std::vector<double> current_state_mean_as_vector{
      data_series.at(data_series.size() - 1) * std::exp(filter_sde.getMu()),
      0.0, 0.0
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

  const BoostMartixInverter matrix_inverter;
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
  for (u_int32_t i{0}; i < current_state_mean.size(); i++)
    posterior_state.current_state_mean(i) = current_state_mean(i);
}
void KcaStates::setCurrentStateMean(std::vector<double>& current_state_mean) {
  move_std_vector_to_vector(
      std::move(current_state_mean), posterior_state.current_state_mean
  );
}
void KcaStates::setCurrentStateCovariance(
    const matrix<double>& current_state_covariance
) {
  for (u_int32_t i{0}; i < current_state_covariance.size1(); i++)
    row(posterior_state.current_state_covariance, i) =
        row(current_state_covariance, i);
}
void KcaStates::setCurrentStateCovariance(
    std::vector<std::vector<double>>& current_state_covariance
) {
  move_std_vectors_to_matrix(
      std::move(current_state_covariance),
      posterior_state.current_state_covariance
  );
}
void KcaStates::setObservationMatrix(const matrix<double>& observation_matrix) {
  for (u_int32_t i{0}; i < observation_matrix.size1(); i++)
    row(prior_state.observation_matrix, i) = row(observation_matrix, i);
}
void KcaStates::setObservationMatrix(
    std::vector<std::vector<double>>& observation_matrix
) {
  move_std_vectors_to_matrix(
      std::move(observation_matrix), prior_state.observation_matrix
  );
}
void KcaStates::setObservationOffset(const double& observation_offset) {
  prior_state.observation_offset = observation_offset;
}
void KcaStates::setPredictedObservationCovariance(
    const matrix<double>& predicted_observation_covariance
) {
  for (u_int32_t i{0}; i < predicted_observation_covariance.size1(); i++)
    row(prior_state.predicted_observation_covariance, i) =
        row(predicted_observation_covariance, i);
}
void KcaStates::setPredictedObservationCovariance(
    std::vector<std::vector<double>>& predicted_observation_covariance
) {
  move_std_vectors_to_matrix(
      std::move(predicted_observation_covariance),
      prior_state.predicted_observation_covariance
  );
}
void KcaStates::setPredictedObservationMean(
    const vector<double>& predicted_observation_mean
) {
  for (u_int32_t i{0}; i < predicted_observation_mean.size(); i++)
    prior_state.predicted_observation_mean(i) = predicted_observation_mean(i);
}
void KcaStates::setPredictedObservationMean(
    std::vector<double>& predicted_observation_mean
) {
  move_std_vector_to_vector(
      std::move(predicted_observation_mean),
      prior_state.predicted_observation_mean
  );
}
void KcaStates::setPredictedStateCovariance(
    const matrix<double>& predicted_state_covariance
) {
  for (u_int32_t i{0}; i < predicted_state_covariance.size1(); i++)
    row(prior_state.predicted_state_covariance, i) =
        row(predicted_state_covariance, i);
}
void KcaStates::setPredictedStateCovariance(
    std::vector<std::vector<double>>& predicted_state_covariance
) {
  move_std_vectors_to_matrix(
      std::move(predicted_state_covariance),
      prior_state.predicted_state_covariance
  );
}
void KcaStates::setPredictedStateMean(
    const vector<double>& predicted_state_mean
) {
  for (u_int32_t i{0}; i < predicted_state_mean.size(); i++)
    prior_state.predicted_state_mean(i) = predicted_state_mean(i);
}
void KcaStates::setPredictedStateMean(
    std::vector<double>& predicted_state_mean
) {
  move_std_vector_to_vector(
      std::move(predicted_state_mean), prior_state.predicted_state_mean
  );
}
void KcaStates::setTransitionCovariance(
    const matrix<double>& transition_covariance
) {
  for (u_int32_t i{0}; i < transition_covariance.size1(); i++)
    row(transition_state.transition_covariance, i) =
        row(transition_covariance, i);
}
void KcaStates::setTransitionCovariance(
    std::vector<std::vector<double>>& transition_covariance
) {
  move_std_vectors_to_matrix(
      std::move(transition_covariance), transition_state.transition_covariance
  );
}
void KcaStates::setTransitionMatrix(const matrix<double>& transition_matrix) {
  for (u_int32_t i{0}; i < transition_matrix.size1(); i++)
    row(transition_state.transition_matrix, i) = row(transition_matrix, i);
}
void KcaStates::setTransitionMatrix(
    std::vector<std::vector<double>>& transition_matrix
) {
  move_std_vectors_to_matrix(
      std::move(transition_matrix), transition_state.transition_matrix
  );
}
