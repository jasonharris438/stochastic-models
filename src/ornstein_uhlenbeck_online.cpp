#include "stochastic_models/likelihood/ornstein_uhlenbeck_online.h"

#include <cmath>

OrnsteinUhlenbeckUpdater::OrnsteinUhlenbeckUpdater(
    const OrnsteinUhlenbeckLikelihoodComponents components,
    const OrnsteinUhlenbeckParameters parameters
)
    : likelihood(OrnsteinUhlenbeckLikelihood()), components(components),
      parameters(parameters) {}

const OrnsteinUhlenbeckParameters
OrnsteinUhlenbeckUpdater::getParameters() const {
  return parameters;
}

const OrnsteinUhlenbeckParameters OrnsteinUhlenbeckUpdater::updateState(
    const double& new_observation, const double& last_observation
) {
  components = likelihood.updateComponents(
      components, new_observation, last_observation
  );
  parameters = likelihood.calculateParameters(components);
  return parameters;
}
