#ifndef STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_LIKELIHOOD_H
#define STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_LIKELIHOOD_H
#include <cstdint>
#include <vector>

/**
 * @brief Maximum-likelihood parameters for the General Linear SDE model.
 *
 * @param mu: estimated autoregressive coefficient (log scale),
 * @param sigma: estimated innovation standard deviation.
 */
struct GeneralLinearParameters {
  double mu;
  double sigma;
};

/**
 * @brief Precomputed summary components used to compute the General Linear SDE
 * model maximum-likelihood estimates.
 *
 * @param lag_squared: sum of squared lagged observations,
 * @param lead_lag_inner_product: inner product of lead and lagged observations,
 * @param squared_error: sum of squared errors between observations and their
 *                       conditional means,
 * @param n_obs: number of observations in the series.
 */
struct GeneralLinearLikelihoodComponents {
  double lag_squared;
  double lead_lag_inner_product;
  double squared_error;
  uint32_t n_obs;
};

class GeneralLinearLikelihoodComponentCalculator {
public:
  /**
   * @brief Safe and convenient division helper to calculate the current mean
   * value of a GeneralLinear SDE model.
   *
   * Calculates the mean value as a component of the mu parameter MLE.
   *
   * @param numerator Numerator value.
   * @param denominator Denominator value (should not be zero in normal use).
   * @return const double Quotient of numerator and denominator.
   */
  const double
  calculateSeriesMean(const double& numerator, const double& denominator) const;
  /**
   * @brief Calculates the inner product of the lead and lag of the data
   * series in `data`.
   *
   * @param data Data series used to calculate.
   * @return const double The inner product of the lead and lag of the data
   * series.
   */
  const double
  calculateLeadLagInnerProduct(const std::vector<double>& data) const;
  /**
   * @brief Calculates the sum of the squares of the lag of the data series in
   * `data`.
   * @param data Data series used to calculate.
   * @return const double The sum of the squares of the lag of the data
   * series.
   */
  const double calculateLagSquared(const std::vector<double>& data) const;
  /**
   * @brief Calculates the sum of squared errors between observations and
   * their conditional means for the data series in `data`.
   *
   * @param data Data series used to calculate.
   * @return const double The sum of squared errors for the data series.
   */
  const double calculateSquaredError(
      const std::vector<double>& data, const double& mu
  ) const;
  /**
   * @brief Calculates the kernel value for the sigma parameter of the data
   * series in `data`.
   * @param n_observations The number of observations in the data series.
   * @param sigma The sigma parameter for the data series.
   * @return const double The kernel value for the sigma parameter.
   */

  /**
   * @brief Incrementally update the lead-lag inner product with one new pair.
   *
   * @param lead_lag_inner_product Current accumulated inner product.
   * @param new_observation New (lead) observation X_t.
   * @param last_observation Previous (lag) observation X_{t-1}.
   * @return const double Updated inner product.
   */
  const double updateLeadLagInnerProduct(
      const double& lead_lag_inner_product,
      const double& new_observation,
      const double& last_observation
  ) const;
  /**
   * @brief Incrementally update the sum of squared lag observations.
   *
   * @param lag_squared Current accumulated lag squared.
   * @param last_observation The new time t-1 observation to add.
   * @return const double Updated lag squared sum.
   */
  const double updateLagSquared(
      const double& lag_squared, const double& last_observation
  ) const;
  /**
   * @brief Update the accumulated squared error term used in sigma MLE.
   *
   * This implements a Welford-like incremental update suitable for
   * one-observation-at-a-time streaming updates.
   *
   * @param squared_error Current accumulated squared error.
   * @param new_observation New time t (lead) observation.
   * @param last_observation Previous time t-1 (lag) observation.
   * @param mean Current running estimate of the mean used in the update.
   * @param n_observations Current number of observations (before this update).
   * @return const double Updated squared error.
   */
  const double updateSquaredError(
      const double& squared_error,
      const double& new_observation,
      const double& last_observation,
      const double& mean,
      const uint32_t& n_observations
  ) const;

  /**
   * @brief Compute the MLE for mu from precomputed components.
   *
   * @param components Precomputed likelihood components for the series.
   * @return const double Estimated mu (log-autoregressive coefficient).
   */
  const double
  calculateMu(const GeneralLinearLikelihoodComponents& components) const;
  /**
   * @brief Compute the MLE for sigma given components and mu.
   *
   * @param components Precomputed likelihood components.
   * @param mu The mu parameter (estimated separately).
   * @return const double Estimated sigma (innovation standard deviation).
   */
  const double calculateSigma(
      const GeneralLinearLikelihoodComponents& components, const double& mu
  ) const;
  /**
   * @brief Calculates the maximum likelihood values for all model
   * parameters.
   *
   * @param data Data series used to calculate.
   * @return const std::unordered_map<std::string, const double> The maximum
   * likelihood value for model parameters.
   */
  /**
   * @brief Calculate the model's conditional variance from its parameters.
   *
   * @param parameters Model parameters (mu, sigma).
   * @return const double Conditional variance of the AR(1)-like model.
   */
  const double
  calculateConditionalVariance(const GeneralLinearParameters& parameters) const;
};

/**
 * @brief Handles calculating maxmimum likelihood parameters specifically under
 * General Linear SDE process assumptions for a data series.
 *
 */
class GeneralLinearLikelihood {
private:
  /**
   * @brief Component calculator instance to use when calculating likelihood
   * components from a data series.
   */
  GeneralLinearLikelihoodComponentCalculator component_calculator;

public:
  /**
   * @brief Calculates the likelihood equation components from data series
   * and model parameters.
   *
   * @param data Data series used to calculate.
   * @return const GeneralLinearLikelihoodComponents The likelihood equation
   * components.
   */
  const GeneralLinearLikelihoodComponents
  calculateComponents(const std::vector<double>& data) const;
  const GeneralLinearLikelihoodComponents updateComponents(
      const GeneralLinearLikelihoodComponents& components,
      const double& new_observation,
      const double& last_observation
  ) const;
  const GeneralLinearParameters calculateParameters(
      const GeneralLinearLikelihoodComponents& components
  ) const;
  const double
  calculateConditionalVariance(const GeneralLinearParameters& parameters) const;
};

#endif // STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_LIKELIHOOD_H
