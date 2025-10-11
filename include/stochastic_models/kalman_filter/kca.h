#ifndef KCA_H
#define KCA_H

#include "stochastic_models/kalman_filter/states.h"

/**
 * @brief Class to expose an interface to use the Kinetic components analysis
 * implementation of the kalman filter to users of this library.
 *
 * Functions are provided to initialise the internal state of the kalman filter
 * with a dataset (vector of values), update the priors and posteriors of the
 * filter and return the current state of the filter.
 *
 * The filter_state member attribute is used to store the current state of the
 * kalman filter. The values provided in instantiation are used to initialise
 * the memory of internally managed attributes.
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
 * @param observation_offset The offset to apply to the observation value.
 */
class KineticComponents {
private:
  KcaStates filter_state;

public:
  KineticComponents(const FilterSystemDimensions& dimensions);
  void setFilterState(const KcaStates& state);
  const KcaStates getFilterState() const;
  const std::vector<double> getCurrentState() const;
  const bool& isInitialised() const;
  const bool& isPriorStateValid() const;
  void initialiseFilter(
      const std::vector<double>& data_series, const double& h, const double& q
  );
  /**
   * @brief Updates the state of the kalman filter prior predicted state mean
   * and covariance.
   *
   * The filter_state member attribute state is updated with new values for
   * the predicted state mean and covariance.
   */
  void updatePriors();
  void
  updatePosteriors(const double& observation, const double& innovation_sigma);

  const std::vector<double> getStateVector() const;
  const std::vector<double> getStandardDevVector() const;
};
#endif // KCA_H
