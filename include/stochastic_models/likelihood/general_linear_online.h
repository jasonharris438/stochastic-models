#ifndef STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_ONLINE_H
#define STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_ONLINE_H

#include "stochastic_models/likelihood/general_linear_likelihood.h"

/**
 * @brief General Linear Model Updater class that handles online updating of
 * the General Linear SDE model specification.
 *
 * This class holds the parameter state required to perform a single-step
 * incremental update of MLE parameters for a general linear SDE model.
 */
class GeneralLinearUpdater {
private:
  /**
   * @brief Stateless likelihood helper instance used for component/parameter
   * calculations.
   */
  const GeneralLinearLikelihood likelihood;
  /**
   * @brief Accumulated parameter components (lag sums, inner-products, squared
   * errors).
   */
  GeneralLinearLikelihoodComponents components;
  /**
   * @brief Current parameter estimates (mu, sigma) used as the baseline for
   * incremental updates.
   */
  GeneralLinearParameters parameters;

public:
  /**
   * @brief Construct updater from precomputed components and parameters.
   *
   * @param components Precomputed components for the observed series.
   * @param parameters Initial parameter estimates (mu, sigma).
   */
  GeneralLinearUpdater(
      const GeneralLinearLikelihoodComponents components,
      const GeneralLinearParameters parameters
  );

  /**
   * @brief Perform a single online update step and return new parameter
   * estimates.
   *
   * This updates the internal `components` and `parameters` using the new
   * observation and returns the updated parameter bundle. The internal state
   * is also updated (the method has side-effects on the stored components).
   *
   * @param new_observation New time t observed value.
   * @param last_observation Previous time t-1 observed value.
   * @return const GeneralLinearParameters Updated parameter bundle.
   */
  const GeneralLinearParameters
  updateState(const double& new_observation, const double& last_observation);
};
#endif // STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_ONLINE_H
