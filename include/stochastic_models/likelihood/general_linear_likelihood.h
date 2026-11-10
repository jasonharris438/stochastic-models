#ifndef STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_LIKELIHOOD_H
#define STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_LIKELIHOOD_H

/**
 * @file include/stochastic_models/likelihood/general_linear_likelihood.h
 * @brief Maximum likelihood estimators for the General Linear SDE model.
 * @ingroup likelihood
 * @details
 * The estimators assume unit spacing between consecutive observations
 * (dt = 1). Applied to a series sampled every t time units of the continuous
 * process they return mu t and sigma sqrt(t), with sigma the diffusion
 * coefficient of the SDE. See docs/derivations/sde-mle-derivations.md.
 */
#include <cstdint>
#include <vector>

/**
 * @brief Maximum-likelihood parameters for the General Linear SDE model.
 *
 * @param mu: estimated autoregressive coefficient (log scale),
 * @param sigma: estimated diffusion coefficient of the SDE.
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
   * @throws InvalidNumberObservationsError If data has fewer than 2 values.
   */
  const double
  calculateLeadLagInnerProduct(const std::vector<double>& data) const;
  /**
   * @brief Calculates the sum of the squares of the lag of the data series in
   * `data`.
   * @param data Data series used to calculate.
   * @return const double The sum of the squares of the lag of the data
   * series.
   * @throws InvalidNumberObservationsError If data has fewer than 2 values.
   */
  const double calculateLagSquared(const std::vector<double>& data) const;
  /**
   * @brief Calculates the sum of squared errors between observations and
   * their conditional means for the data series in `data`.
   *
   * @param data Data series used to calculate.
   * @return const double The sum of squared errors for the data series.
   * @throws InvalidNumberObservationsError If data has fewer than 2 values.
   */
  const double calculateSquaredError(
      const std::vector<double>& data, const double& mu
  ) const;
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
   * @brief Adds 1 pair to the accumulated squared error exactly.
   *
   * Adds `(new - prior_mean * last)^2 * lag_squared_old / lag_squared_new`,
   * which equals a batch recalculation over the extended series. A zero old
   * sum means no prior pair, and adds 0.
   *
   * @param squared_error Current accumulated squared error.
   * @param new_observation New time t (lead) observation.
   * @param last_observation Previous time t-1 (lag) observation.
   * @param prior_mean The autoregressive coefficient before this pair.
   * @param lag_squared_old The sum of squared lags before this pair.
   * @param lag_squared_new The sum of squared lags including this pair.
   * @return const double Updated squared error.
   */
  const double updateSquaredError(
      const double& squared_error,
      const double& new_observation,
      const double& last_observation,
      const double& prior_mean,
      const double& lag_squared_old,
      const double& lag_squared_new
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
   * @return const double Estimated sigma, the diffusion coefficient of the
   * SDE, from the squared error over the n - 1 pairs.
   */
  const double calculateSigma(
      const GeneralLinearLikelihoodComponents& components, const double& mu
  ) const;
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
   * @throws InvalidNumberObservationsError If `data.size() < 2` (at least one
   *         lead/lag pair is required).
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
};

#endif // STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_LIKELIHOOD_H
