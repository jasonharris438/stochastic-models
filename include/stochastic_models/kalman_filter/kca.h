#ifndef STOCHASTIC_MODELS_KALMAN_FILTER_KCA_H
#define STOCHASTIC_MODELS_KALMAN_FILTER_KCA_H

#include "stochastic_models/kalman_filter/states.h"

/**
 * @file
 * @brief High-level Kinetic Components Analysis (KCA) Kalman filter facade.
 */

/**
 * @brief Facade providing a simple interface to initialize and run the
 * KCA-style Kalman filter.
 *
 * The class wraps the lower-level `KcaStates` and provides methods that perform
 * functions such as: initialize, update priors, update posteriors and accessors
 * for the current state.
 */
class KineticComponents {
private:
  KcaStates filter_state;

public:
  /**
   * @brief Construct a KineticComponents instance with the given dimensions.
   * @param dimensions Dimensions that determine internal buffer sizes.
   */
  KineticComponents(const FilterSystemDimensions& dimensions);

  /**
   * @brief Replace the internal `KcaStates` with an externally-constructed
   * instance.
   */
  void setFilterState(const KcaStates& state);

  /**
   * @brief Return a copy of the internal filter state.
   */
  const KcaStates getFilterState() const;

  /**
   * @brief Return the current state mean as a std::vector (copy).
   */
  const std::vector<double> getCurrentState() const;

  /**
   * @brief Query whether the internal filter has been initialised.
   */
  const bool& isInitialised() const;

  /**
   * @brief Query whether priors have been successfully computed and are valid.
   */
  const bool& isPriorStateValid() const;

  /**
   * @brief Initialize the filter using a raw data series and numeric params.
   *
   * @param data_series Series used to estimate initial state.
   * @param h finite-difference spacing used by the KCA estimator.
   * @param q transition covariance scaling.
   */
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
#endif // STOCHASTIC_MODELS_KALMAN_FILTER_KCA_H
