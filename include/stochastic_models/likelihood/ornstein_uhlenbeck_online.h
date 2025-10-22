#ifndef STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_ONLINE_H
#define STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_ONLINE_H

#include "stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h"

/**
 * @file
 * @brief Stateful updater for Ornstein-Uhlenbeck MLE parameters.
 */

/**
 * @brief Online updater for Ornstein-Uhlenbeck parameters.
 *
 * Holds a copy of the current component summary and parameter bundle and
 * exposes an update step that accepts a new observation pair.
 */
class OrnsteinUhlenbeckUpdater {
private:
  const OrnsteinUhlenbeckLikelihood likelihood;
  OrnsteinUhlenbeckLikelihoodComponents components;
  OrnsteinUhlenbeckParameters parameters;

public:
  /**
   * @brief Construct updater from components and parameters (copied).
   */
  OrnsteinUhlenbeckUpdater(
      const OrnsteinUhlenbeckLikelihoodComponents components,
      const OrnsteinUhlenbeckParameters parameters
  );

  /**
   * @brief Return the current parameter bundle.
   */
  const OrnsteinUhlenbeckParameters getParameters() const;

  /**
   * @brief Update internal components using a new observation and return
   * updated parameters.
   */
  const OrnsteinUhlenbeckParameters
  updateState(const double& new_observation, const double& last_observation);
};

#endif // STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_ONLINE_H
