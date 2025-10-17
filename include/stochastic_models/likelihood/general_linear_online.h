#ifndef STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_ONLINE_H
#define STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_ONLINE_H

/**
 * @brief General Linear Model Updater class that handles online updating of
 * the General Linear SDE model specification.
 *
 * Welford's (1962) method for calculating corrected sums of squares is used to
 * update the model parameters.
 */
class GeneralLinearUpdater {
private:
  /**
   * @brief The numerator used when tracking changes to the mu model parameter
   * internally.
   */
  double mu_numerator;
  /**
   * @brief The denominator used when tracking changes to the mu model parameter
   * internally.
   */
  double mu_denominator;
  /**
   * @brief The state value used when tracking changes to the sigma model
   * parameter internally.
   */
  double sigma_state_value;
  /**
   * @brief The last observation used to perform a model update.
   */
  double last_observation;
  /**
   * @brief The last number of observations used to perform a model update.
   */
  int n_obs;

  /**
   * @brief Update internal values tracked for online calculation of the
   * expected value parameter state.
   * @param mu_numerator_new The new numerator value used in the expected value
   * calculation.
   * @param mu_denominator_new The new denominator value used in the expected
   * value calculation.
   */
  void applyInternalMuUpdate(
      const double mu_numerator_new, const double mu_denominator_new
  );

  /**
   * @brief Update internal values tracked for online calculation of the
   * standard deviation parameter state.
   * @param sigma_state_value_new The new state value used in the standard
   * deviation calculation.
   */
  void applyInternalSigmaUpdate(const double sigma_state_value_new);

public:
  GeneralLinearUpdater(
      const double mu_numerator,
      const double mu_denominator,
      const double sigma_state_value,
      const double last_observation,
      const int n_obs
  );

  /**
   * @brief Returns the expected value of the observation given the current
   * model state.
   * @param observation The x_t-1 observation used to calculate the expected
   * value.
   * @return The expected value of the observation given the current model
   * state.
   */
  const double getExpectedValue(const double& observation) const;

  /**
   * @brief Returns the current value of the expected value model parameter.
   * @return The current value of the expected value model parameter.
   */
  const double getMu() const;

  /**
   * @brief Returns the current value of the standard deviation model parameter.
   * @return The current value of the standard deviation model parameter.
   */
  const double getSigma() const;

  /**
   * @brief Returns the current number of observations used in the model.
   * @return The current number of observations used in the model.
   */
  const int getNumberObs() const;

  /**
   * @brief Update the expected value model parameter using a new observation.
   * @param value The new observation used to update the expected value model
   * parameter.
   * @return The updated expected value model parameter.
   */
  const double updateMu(const double& value);

  /**
   * @brief Update the standard deviation model parameter using a new
   * observation.
   * @param value The new observation used to update the standard deviation
   * model parameter.
   * @return The updated standard deviation model parameter.
   */
  const double updateSigma(const double& value);

  /**
   * @brief Update internal values tracked for online calculation of the general
   * internal state of the model.
   * @param n_obs_new The new number of observations used in the general
   * internal state calculation.
   * @param new_observation The new observation used in the general internal
   * state calculation.
   */
  void internalStateUpdate(const int n_obs_new, const double new_observation);
};
#endif // STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_ONLINE_H
