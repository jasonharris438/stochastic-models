#ifndef STATES_H
#define STATES_H
#include "stochastic_models/likelihood/general_linear_likelihood.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

// Just for this module as we do not introduce any other namespaces.
using namespace boost::numeric::ublas;

/**
 * @brief Class to represent the predicted state internals of the Kalman Filter.
 *
 * This class is used to calculate the prior predicted state of the Kalman
 * Filter before a new observation is received.
 *
 * Methods are implemented to calculate and return the mean and covariance of
 * the predicted state.
 *
 * @param transition_matrix The transition matrix of the Kalman Filter.
 *
 */
class PredictedState {
private:
  matrix<double> transition_matrix;

public:
  /**
   * @brief Constructor for the PredictedState class.
   *
   * This constructor accepts the transition matrix of the predicted
   * state as a boost uBLAS object constructed externally.
   *
   * @param transition_matrix The transition matrix of the Kalman Filter.
   */
  PredictedState(matrix<double> transition_matrix);
  /**
   * @brief Method to set the transition matrix of the prior predicted state.
   * @param transition_matrix The transition matrix for the predicted state.
   */
  void setTransitionMatrix(const matrix<double>& transition_matrix);
  /**
   * @brief Method to calculate the mean of the predicted state.
   * @param current_state_mean The mean of the current state.
   * @return The mean of the predicted state.
   */
  const vector<double>
  calculateMean(const vector<double>& current_state_mean) const;
  /**
   * @brief Method to calculate the covariance of the predicted state.
   * @param current_state_covariance The covariance of the current state.
   * @param transition_covariance The covariance of the transition matrix.
   * @return The covariance of the predicted state.
   */
  const matrix<double> calculateCovariance(
      const matrix<double>& current_state_covariance,
      const matrix<double>& transition_covariance
  ) const;
};

/**
 * @brief Class to represent the prior predicted observation internals of
 * the Kalman Filter.
 *
 * This class is used to calculate the predicted observation mean and covariance
 * for a Kalman Filter to be made before a new observation is received.
 *
 * @param observation_matrix The predicted observation matrix of the Kalman
 * Filter.
 * @param observation_offset The offset value to add to the predicted
 * observation mean.
 *
 */
class PredictedObservation {
private:
  matrix<double> observation_matrix;
  double observation_offset;

public:
  /**
   * @brief Constructor for the PredictedObservation class.
   *
   * This constructor accepts the observation matrix and offset of the
   * predicted observation as uBLAS objects constructed externally.
   *
   */
  PredictedObservation(
      matrix<double> observation_matrix, double observation_offset
  );
  /**
   * @brief Method to set the observation matrix of the prior predicted
   * observation.
   * @param observation_matrix The observation matrix of the predicted
   * observation.
   */
  void setObservationMatrix(const matrix<double>& observation_matrix);
  const matrix<double>& getObservationMatrix() const;
  /**
   * @brief Method to set the observation offset of the prior predicted
   * observation.
   * @param observation_offset The offset value to add to the predicted
   * observation mean.
   */
  void setObservationOffset(const double& observation_offset);
  /**
   * @brief Method to calculate the mean of the prior predicted observation.
   * @param predicted_state_mean The mean of the prior predicted state.
   * @return The mean of the prior predicted observation.
   */
  const vector<double>
  calculateMean(const vector<double>& predicted_state_mean) const;
  /**
   * @brief Method to calculate the covariance of the prior predicted
   * observation.
   * @param predicted_state_covariance The covariance of the prior predicted
   * state.
   * @param innovation_sigma The sigma parameter of the observation
   * innovation.
   * @return The covariance of the prior predicted observation.
   */
  const matrix<double> calculateCovariance(
      const matrix<double>& predicted_state_covariance,
      const double& innovation_sigma
  ) const;
  /**
   * @brief Method to calculate the kalman gain of the predicted state.
   * @param predicted_state_covariance The covariance of the prior predicted
   * state.
   * @param predicted_observation_covariance The covariance of the prior
   * predicted observation.
   * @return The kalman gain of the predicted state.
   */
  template <typename T>
  const matrix<double> calculateKalmanGain(
      const matrix<double>& predicted_state_covariance,
      const matrix<double>& predicted_observation_covariance,
      const T& matrix_inverter
  ) const {
    const matrix<double> observation_transpose = trans(observation_matrix);
    const matrix<double> inverse_covariance =
        matrix_inverter.invertMatrix(predicted_observation_covariance);
    const matrix<double> inner_matrix_product =
        prod(observation_transpose, inverse_covariance);
    const matrix<double> result =
        prod(predicted_state_covariance, inner_matrix_product);
    return result;
  };
};

/**
 * @brief Class to represent the current state internals of the Kalman Filter.
 *
 * This class is used to calculate the predicted posterior current state of the
 * Kalman Filter after a new observation is received.
 *
 * Methods are implemented to calculate and return the mean and covariance of
 * the posterior current state.
 *
 */
class CurrentState {
public:
  /**
   * @brief Method to calculate the predicted mean of the posterior current
   * state.
   * @param predicted_state_mean The mean of the prior predicted state.
   * @param kalman_gain The Kalman gain of the Kalman Filter.
   * @param innovation The innovation (difference between predicted and actual
   * observation values) of the Kalman Filter.
   * @return The predicted mean of the posterior current state.
   */
  const vector<double> calculateMean(
      const vector<double>& predicted_state_mean,
      const vector<double>& kalman_gain,
      const double& innovation
  ) const;
  /**
   * @brief Method to calculate the predicted covariance of the posterior
   * current state.
   * @param predicted_state_covariance The covariance of the prior predicted
   * state.
   * @param observation_matrix The observation matrix of the Kalman Filter.
   * @param kalman_gain The last calculated Kalman gain of the Kalman Filter
   * provided as a matrix.
   */
  const matrix<double> calculateCovariance(
      const matrix<double>& predicted_state_covariance,
      const matrix<double>& observation_matrix,
      const matrix<double>& kalman_gain
  ) const;
};

/**
 * @brief Struct to represent the prior state of the Kalman Filter.
 *
 * Contains the predicted observation and state mean and covariance matrices.
 *
 * @param state_mean_dimension The dimension of the state mean vector.
 * @param state_covariance_rows The number of rows in the state covariance
 * matrix.
 * @param state_covariance_columns The number of columns in the state covariance
 * matrix.
 * @param observation_matrix_rows The number of rows in the observation matrix.
 * @param observation_matrix_columns The number of columns in the observation
 * matrix.
 * @param observation_covariance_rows The number of rows in the observation
 * covariance matrix.
 * @param observation_covariance_columns The number of columns in the
 * observation covariance matrix.
 * @param observation_offset The offset value to add to the predicted
 * observation mean.
 *
 */
struct PriorState {
  vector<double> predicted_observation_mean;
  vector<double> predicted_state_mean;
  matrix<double> predicted_observation_covariance;
  matrix<double> predicted_state_covariance;
  matrix<double> observation_matrix;
  double observation_offset;

  PriorState(
      const int& state_mean_dimension,
      const int& state_covariance_rows,
      const int& state_covariance_columns,
      const int& observation_matrix_rows,
      const int& observation_matrix_columns,
      const int& observation_covariance_rows,
      const int& observation_covariance_columns,
      const double& observation_offset
  );
};

/**
 * @brief Struct to represent the posterior state of the Kalman Filter.
 *
 * Contains the current state mean and covariance matrices.
 *
 * @param current_state_mean The current state mean vector.
 * @param current_state_covariance The current state covariance matrix.
 */
struct PosteriorState {
  vector<double> current_state_mean;
  matrix<double> current_state_covariance;

  PosteriorState(
      const int& state_mean_dimension,
      const int& state_covariance_rows,
      const int& state_covariance_columns
  );
};

/**
 * @brief Struct to represent the transition state of the Kalman Filter.
 *
 * Contains components that determine the transition between states. The
 * transition matrix and covariance matrices are stored.
 *
 * @param transition_matrix The transition matrix of the Kalman Filter.
 * @param transition_covariance The transition covariance matrix of the Kalman
 * Filter.
 */
struct TransitionState {
  matrix<double> transition_matrix;
  matrix<double> transition_covariance;

  TransitionState(
      const int& state_covariance_rows, const int& state_covariance_columns
  );
};

/**
 * @brief Struct to represent the boolean state of the Kalman Filter.
 *
 * Contains boolean flags to determine the state of the Kalman Filter.
 *
 * @param initialised Flag to determine if the Kalman Filter has been
 * initialised.
 * @param priors_set Flag to determine if the Kalman Filter priors have been
 * set.
 */
struct FilterState {
  bool initialised;
  bool priors_set;

  FilterState();
};

/**
 * @brief Handles internals of the stochastic differential equation (SDE) that
 * governs the process being analysed by the kalman filtering process within the
 * KCA.
 *
 * The maximum likelihood parameters are calculated and stored within this
 * class.
 *
 * @param likelihood The likelihood calculator for the SDE process.
 * @param mu The mu parameter of the SDE process.
 * @param sigma The sigma parameter of the SDE process.
 * @param conditional_variance The conditional variance of the SDE process.
 * @param mu_numerator An internally managed value used to calculate the mu
 * parameter in initialisation and on update.
 * @param mu_denominator An internally managed value used to calculate the mu
 * parameter in initialisation and on update.
 * @param sigma_kernel_value An internally managed value used to calculate the
 * sigma parameter in initialisation and on update.
 * @param n_observations The number of observations in the data series used to
 * calculate the likelihood parameters.
 */
class FilterGeneralSde {
private:
  GeneralLinearLikelihood likelihood;
  double mu;
  double sigma;
  double conditional_variance;
  double mu_numerator;
  double mu_denominator;
  double sigma_kernel_value;
  u_int32_t n_observations;

public:
  FilterGeneralSde();
  FilterGeneralSde(
      const double& mu,
      const double& sigma,
      const double& conditional_variance,
      const double& mu_numerator,
      const double& mu_denominator,
      const double& sigma_kernel_value,
      const u_int32_t& n_observations
  );
  /**
   * @brief Retrieves the current value of the mu parameter from the SDE
   * process.
   * @return The mu parameter of the SDE process.
   */
  const double& getMu() const;
  /**
   * @brief Retrieves the current value of the sigma parameter from the SDE
   * process.
   * @return The sigma parameter of the SDE process.
   */
  const double& getSigma() const;
  /**
   * @brief Retrieves the current value of the conditional covariance from the
   * SDE process.
   * @return The conditional covariance parameter of the SDE process.
   */
  const double& getConditionalVariance() const;
  /**
   * @brief Initializes completely the likelihood state of the SDE process
   * from a data series.
   * @param data_series The data series used to calculate the likelihood
   * parameters.
   */
  void initializeLikelihoodState(const std::vector<double>& data_series);
};

/**
 * @brief Contains the dimensions of a Kalman Filter system.
 *
 * Internal system state dimensions are stored in this struct.
 *
 * @param state_mean_dimension The dimension of the state mean vector.
 * @param state_covariance_rows The number of rows in the state covariance
 * matrix.
 * @param state_covariance_columns The number of columns in the state covariance
 * matrix.
 * @param observation_matrix_rows The number of rows in the observation matrix.
 * @param observation_matrix_columns The number of columns in the observation
 * matrix.
 * @param observation_covariance_rows The number of rows in the observation
 * covariance matrix.
 * @param observation_covariance_columns The number of columns in the
 * observation covariance matrix.
 * @param observation_offset The offset value to add to the predicted
 * observation mean.
 */
struct FilterSystemDimensions {
  int state_mean_dimension;
  int state_covariance_rows;
  int state_covariance_columns;
  int observation_matrix_rows;
  int observation_matrix_columns;
  int observation_covariance_rows;
  int observation_covariance_columns;
  double observation_offset;

  FilterSystemDimensions();
  FilterSystemDimensions(
      int state_mean_dimension,
      int state_covariance_rows,
      int state_covariance_columns,
      int observation_matrix_rows,
      int observation_matrix_columns,
      int observation_covariance_rows,
      int observation_covariance_columns,
      double observation_offset
  );
};

/**
 * @brief Class to represent the state handler for the kinetic components
 * analysis (KCA) implementation.
 *
 * The complete state is stored and managed by this class, separating it from
 * the low-level detail of the kalman filter.
 *
 * @param dimensions The dimensions of the internals of the Kalman Filter
 * system.
 */
class KcaStates {
private:
  PriorState prior_state;
  PosteriorState posterior_state;
  TransitionState transition_state;
  FilterState filter_state;
  FilterGeneralSde filter_sde;

  /**
   * @brief Moves a std::vector of std::vectors to a boost uBLAS matrix.
   * @param matrix_as_vectors The matrix represented as a vector of vectors.
   * @param target The target boost matrix to move the vectors to.
   */
  void move_std_vectors_to_matrix(
      std::vector<std::vector<double>>&& matrix_as_vectors,
      matrix<double>& target
  );
  /**
   * @brief Moves a std::vector to a boost uBLAS vector.
   * @param vector_as_vector The vector represented as a vector.
   * @param target The target boost vector to move the vector to.
   */
  void move_std_vector_to_vector(
      std::vector<double>&& vector_as_vector, vector<double>& target
  );

public:
  KcaStates(const FilterSystemDimensions& dimensions);

  /**
   * @brief Initializes the initial KCA system state by setting the starting
   * values of the internally managed state objects with data derived from the
   * data_series value provided.
   * @param data_series A vector containing a series of values to base the
   * initialization of the system internals on.
   * @param h A value determining the first and second derivative values of
   * the KCA system.
   * @param q A value determining the transition covariance of the KCA system.
   */
  void setInitialState(
      const std::vector<double>& data_series, const double& h, const double& q
  );
  /**
   * @brief Makes the prediction step and updates the prior predicted state of
   * the KCA system.
   *
   * The system must be initialized before this method is called.
   */
  void updatePredictedState();
  /**
   * @brief Makes the current state update step given observed data and
   * updates the posterior current state of the KCA system.
   * @param observation The observed data value.
   * @param innovation_sigma The sigma value of the innovation of the observed
   * data.
   */
  void
  updateCurrentState(const double& observation, const double& innovation_sigma);
  /**
   * @brief Retrieves the current state mean vector from the KCA system as a
   * std::vector.
   * @return The current state mean vector of the KCA system as a std::vector.
   */
  const std::vector<double> getCurrentStateMeanVector() const;
  /**
   * @brief Retrieves the current state mean from the KCA system.
   * @return The current state mean of the KCA system.
   */
  const vector<double>& getCurrentStateMean() const;
  /**
   * @brief Retrieves the current state covariance from the KCA system.
   * @return The current state covariance of the KCA system.
   */
  const matrix<double>& getCurrentStateCovariance() const;
  /**
   * @brief Retrieves the observation matrix from the KCA system.
   * @return The observation matrix of the KCA system.
   */
  const matrix<double>& getObservationMatrix() const;
  /**
   * @brief Retrieves the observation offset from the KCA system.
   * @return The observation offset of the KCA system.
   */
  const double& getObservationOffset() const;
  /**
   * @brief Retrieves the predicted observation covariance from the KCA
   * system.
   * @return The predicted observation covariance of the KCA system.
   */
  const matrix<double>& getPredictedObservationCovariance() const;
  /**
   * @brief Retrieves the predicted observation mean from the KCA system.
   * @return The predicted observation mean of the KCA system.
   */
  const vector<double>& getPredictedObservationMean() const;
  /**
   * @brief Retrieves the predicted state covariance from the KCA system.
   * @return The predicted state covariance of the KCA system.
   */
  const matrix<double>& getPredictedStateCovariance() const;
  /**
   * @brief Retrieves the predicted state mean from the KCA system.
   * @return The predicted state mean of the KCA system.
   */
  const vector<double>& getPredictedStateMean() const;
  /**
   * @brief Retrieves the transition covariance from the KCA system.
   * @return The transition covariance of the KCA system.
   */
  const matrix<double>& getTransitionCovariance() const;
  /**
   * @brief Retrieves the transition matrix from the KCA system.
   * @return The transition matrix of the KCA system.
   */
  const matrix<double>& getTransitionMatrix() const;

  /**
   * @brief Retrieves a flag indicating the current initialization state of
   * the KCA system.
   * @return A boolean flag indicating whether the system has been initialized
   * and is valid.
   */
  const bool& isInitialised() const;
  /**
   * @brief Retrieves a flag indicating the current state of the priors of the
   * KCA system.
   *
   * If the flag is true, then the posterior state can be updated.
   * @return A boolean flag indicating whether the priors of the system are
   * in a valid state.
   */
  const bool& arePriorsValid() const;
  /**
   * @brief Sets the 'is initialized' flag of the KCA system to true.
   */
  void setInitialized();
  /**
   * @brief Sets the 'are priors valid' flag of the KCA system to true.
   */
  void setPriorsTrue();
  /**
   * @brief Sets the 'are priors valid' flag of the KCA system to false.
   */
  void setPriorsFalse();
  /**
   * @brief Sets the current state mean of the KCA system with a boost uBLAS
   * vector.
   *
   * The values of the vector are copied to the system's current state mean
   * object.
   * @param current_state_mean The current state mean uBLAS vector to copy
   * into the system's current state mean.
   */
  void setCurrentStateMean(const vector<double>& current_state_mean);
  /**
   * @brief Sets the current state mean of the KCA system with a std::vector.
   *
   * The values of the vector are moved to the system's current state mean
   * object.
   * @param current_state_mean The current state mean std::vector to move into
   * the system's current state mean.
   */
  void setCurrentStateMean(std::vector<double>& current_state_mean);
  /**
   * @brief Sets the current state covariance of the KCA system with a boost
   * uBLAS matrix.
   *
   * The values of the matrix are copied to the system's current state
   * covariance object.
   * @param current_state_covariance The current state covariance uBLAS matrix
   * to copy into the system's current state covariance.
   */
  void
  setCurrentStateCovariance(const matrix<double>& current_state_covariance);
  /**
   * @brief Sets the current state covariance of the KCA system with a
   * std::vector of std::vectors.
   *
   * The values of the vector are moved to the system's current state
   * covariance object.
   * @param current_state_covariance The current state covariance std::vector
   * of std::vectors to move into the system's current state covariance.
   */
  void setCurrentStateCovariance(
      std::vector<std::vector<double>>& current_state_covariance
  );
  /**
   * @brief Sets the observation matrix of the KCA system with a boost uBLAS
   * matrix.
   *
   * The values of the matrix are copied to the system's observation matrix
   * object.
   * @param observation_matrix The observation matrix uBLAS matrix to copy
   * into the system's observation matrix.
   */
  void setObservationMatrix(const matrix<double>& observation_matrix);
  /**
   * @brief Sets the observation matrix of the KCA system with a std::vector
   * of std::vectors.
   *
   * The values of the vector are moved to the system's observation matrix
   * object.
   * @param observation_matrix The observation matrix std::vector of
   * std::vectors to move into the system's observation matrix.
   */
  void
  setObservationMatrix(std::vector<std::vector<double>>& observation_matrix);
  /**
   * @brief Sets the observation offset value of the KCA system.
   * @param observation_offset The observation offset value to set.
   */
  void setObservationOffset(const double& observation_offset);
  /**
   * @brief Sets the predicted observation covariance of the KCA system with a
   * boost uBLAS matrix.
   *
   * The values of the matrix are copied to the system's predicted observation
   * covariance object.
   * @param predicted_observation_covariance The predicted observation
   * covariance uBLAS matrix to copy into the system's predicted observation
   * covariance.
   */
  void setPredictedObservationCovariance(
      const matrix<double>& predicted_observation_covariance
  );
  /**
   * @brief Sets the predicted observation covariance of the KCA system with a
   * std::vector of std::vectors.
   *
   * The values of the vector are moved to the system's predicted observation
   * covariance object.
   * @param predicted_observation_covariance The predicted observation
   * covariance std::vector of std::vectors to move into the system's
   * predicted observation covariance.
   */
  void setPredictedObservationCovariance(
      std::vector<std::vector<double>>& predicted_observation_covariance
  );
  /**
   * @brief Sets the predicted observation mean of the KCA system with a boost
   * uBLAS vector.
   *
   * The values of the vector are copied to the system's predicted observation
   * mean object.
   * @param predicted_observation_mean The predicted observation mean uBLAS
   * vector to copy into the system's predicted observation mean.
   */
  void
  setPredictedObservationMean(const vector<double>& predicted_observation_mean);
  /**
   * @brief Sets the predicted observation mean of the KCA system with a
   * std::vector.
   *
   * The values of the vector are moved to the system's predicted observation
   * mean object.
   * @param predicted_observation_mean The predicted observation mean
   * std::vector to move into the system's predicted observation mean.
   */
  void
  setPredictedObservationMean(std::vector<double>& predicted_observation_mean);
  /**
   * @brief Sets the predicted state covariance of the KCA system with a boost
   * uBLAS matrix.
   *
   * The values of the matrix are copied to the system's predicted state
   * covariance object.
   * @param predicted_state_covariance The predicted state covariance uBLAS
   * matrix to copy into the system's predicted state covariance.
   */
  void
  setPredictedStateCovariance(const matrix<double>& predicted_state_covariance);
  /**
   * @brief Sets the predicted state covariance of the KCA system with a
   * std::vector of std::vectors.
   *
   * The values of the vector are moved to the system's predicted state
   * covariance object.
   * @param predicted_state_covariance The predicted state covariance
   * std::vector of std::vectors to move into the system's predicted state
   * covariance.
   */
  void setPredictedStateCovariance(
      std::vector<std::vector<double>>& predicted_state_covariance
  );
  /**
   * @brief Sets the predicted state mean of the KCA system with a boost uBLAS
   * vector.
   *
   * The values of the vector are copied to the system's predicted state mean
   * object.
   * @param predicted_state_mean The predicted state mean uBLAS vector to copy
   * into the system's predicted state mean.
   */
  void setPredictedStateMean(const vector<double>& predicted_state_mean);
  /**
   * @brief Sets the predicted state mean of the KCA system with a
   * std::vector.
   *
   * The values of the vector are moved to the system's predicted state mean
   * object.
   * @param predicted_state_mean The predicted state mean std::vector to move
   * into the system's predicted state mean.
   */
  void setPredictedStateMean(std::vector<double>& predicted_state_mean);
  /**
   * @brief Sets the transition covariance of the KCA system with a boost
   * uBLAS matrix.
   *
   * The values of the matrix are copied to the system's transition covariance
   * object.
   * @param transition_covariance The transition covariance uBLAS matrix to
   * copy into the system's transition covariance.
   */
  void setTransitionCovariance(const matrix<double>& transition_covariance);
  /**
   * @brief Sets the transition covariance of the KCA system with a
   * std::vector of std::vectors.
   *
   * The values of the vector are moved to the system's transition covariance
   * object.
   * @param transition_covariance The transition covariance std::vector of
   * std::vectors to move into the system's transition covariance.
   */
  void setTransitionCovariance(
      std::vector<std::vector<double>>& transition_covariance
  );
  /**
   * @brief Sets the transition matrix of the KCA system with a boost uBLAS
   * matrix.
   *
   * The values of the matrix are copied to the system's transition matrix
   * object.
   * @param transition_matrix The transition matrix uBLAS matrix to copy into
   * the system's transition matrix.
   */
  void setTransitionMatrix(const matrix<double>& transition_matrix);
  /**
   * @brief Sets the transition matrix of the KCA system with a std::vector of
   * std::vectors.
   *
   * The values of the vector are moved to the system's transition matrix
   * object.
   * @param transition_matrix The transition matrix std::vector of
   * std::vectors to move into the system's transition matrix.
   */
  void setTransitionMatrix(std::vector<std::vector<double>>& transition_matrix);
};
#endif // STATES_H
