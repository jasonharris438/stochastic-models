#ifndef STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_LIKELIHOOD_H
#define STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_LIKELIHOOD_H
#include <cstdint>
#include <vector>

struct GeneralLinearParameters {
  double mu;
  double sigma;
};

struct GeneralLinearLikelihoodComponents {
  double lag_squared;
  double lead_lag_inner_product;
  double squared_error;
  uint32_t n_obs;
};

class GeneralLinearLikelihoodComponentCalculator {
public:
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
  const double calculateSigmaKernel(
      const GeneralLinearLikelihoodComponents& components,
      const GeneralLinearParameters& parameters
  ) const;

  const double updateLeadLagInnerProduct(
      const double& lead_lag_inner_product,
      const double& new_observation,
      const double& last_observation
  ) const;
  const double updateLagSquared(
      const double& lag_squared, const double& last_observation
  ) const;
  const double updateSquaredError(
      const double& squared_error,
      const double& new_observation,
      const double& last_observation,
      const double& mean,
      const uint32_t& n_observations
  ) const;

  /**
   * @brief Calculates the `mu` parameter for the data series in `data`.
   * @return const double The `mu` parameter for the data series.
   */
  const double
  calculateMu(const GeneralLinearLikelihoodComponents& components) const;
  /**
   * @brief Calculates the `sigma` parameter for the data series in `data`.
   * @param mu The `mu` parameter for the data series.
   * @return const double The `sigma` parameter for the data series.
   */
  const double calculateSigma(
      const GeneralLinearLikelihoodComponents& components, const double& mu
  ) const;
  /**
   * @brief Calculates the maxmimum likelihood values for all model
   * parameters.
   *
   * @param data Data series used to calculate.
   * @return const std::unordered_map<std::string, const double> The maximum
   * likelihood value for model parameters.
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
  const double calculateSigmaKernel(
      const GeneralLinearLikelihoodComponents& components,
      const GeneralLinearParameters& parameters
  ) const;
};

#endif // STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_LIKELIHOOD_H
