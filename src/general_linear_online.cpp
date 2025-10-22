#include "stochastic_models/likelihood/general_linear_online.h"

#include <cmath>
GeneralLinearUpdater::GeneralLinearUpdater(
    const GeneralLinearLikelihoodComponents components,
    const GeneralLinearParameters parameters
)
    : components(components), parameters(parameters),
      likelihood(GeneralLinearLikelihood()) {}

const GeneralLinearParameters GeneralLinearUpdater::updateState(
    const double& new_observation, const double& last_observation
) {
  components = likelihood.updateComponents(
      components, new_observation, last_observation
  );
  parameters = likelihood.calculateParameters(components);
  return parameters;
}
