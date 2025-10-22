#ifndef STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_ONLINE_H
#define STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_ONLINE_H

#include "stochastic_models/likelihood/general_linear_likelihood.h"

/**
 * @brief General Linear Model Updater class that handles online updating of
 * the General Linear SDE model specification.
 *
 * Welford's (1962) method for calculating corrected sums of squares is used to
 * update the model parameters.
 */
class GeneralLinearUpdater {
private:
  const GeneralLinearLikelihood likelihood;
  GeneralLinearLikelihoodComponents components;
  GeneralLinearParameters parameters;

public:
  GeneralLinearUpdater(
      const GeneralLinearLikelihoodComponents components,
      const GeneralLinearParameters parameters
  );

  /**
   * @brief Update internal values tracked for online calculation of the general
   * internal state of the model.
   * @param new_observation The new observation used in the general internal
   * state calculation.
   */
  const GeneralLinearParameters
  updateState(const double& new_observation, const double& last_observation);
};
#endif // STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_ONLINE_H
