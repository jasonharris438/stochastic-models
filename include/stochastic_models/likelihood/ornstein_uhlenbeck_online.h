#ifndef STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_ONLINE_H
#define STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_ONLINE_H

#include "stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h"

class OrnsteinUhlenbeckUpdater {
private:
  const OrnsteinUhlenbeckLikelihood likelihood;
  OrnsteinUhlenbeckLikelihoodComponents components;
  OrnsteinUhlenbeckParameters parameters;

public:
  OrnsteinUhlenbeckUpdater(
      const OrnsteinUhlenbeckLikelihoodComponents components,
      const OrnsteinUhlenbeckParameters parameters
  );

  const OrnsteinUhlenbeckParameters getParameters() const;
  const OrnsteinUhlenbeckParameters
  updateState(const double& new_observation, const double& last_observation);
};

#endif // STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_ONLINE_H
