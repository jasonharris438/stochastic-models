#ifndef STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_LIKELIHOOD_H
#define STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_LIKELIHOOD_H
#include <cstdint>
#include <vector>

struct GeneralLinearParameters {
  double mu;
  double sigma;
};

class GeneralLinearLikelihoodComponents {
public:
  const double lag_squared;
  const double lead_lag_inner_product;
  const double sigma_kernel;
  const uint32_t n_obs;

  GeneralLinearLikelihoodComponents(
      const double lag_squared,
      const double lead_lag_inner_product,
      const double sigma_kernel,
      const uint32_t n_obs
  );
};

/**
 * @brief Handles calculating maxmimum likelihood parameters specifically under
 * General Linear SDE process assumptions for a data series.
 *
 */
class GeneralLinearLikelihood {
private:
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

public:
  /**
   * @brief Calculates the kernel value for the sigma parameter of the data
   * series in `data`.
   * @param n_observations The number of observations in the data series.
   * @param sigma The sigma parameter for the data series.
   * @return const double The kernel value for the sigma parameter.
   */
  const double calculateSigmaKernel(
      const uint32_t& n_observations, const double& sigma
  ) const;

  /**
   * @brief Calculates the `mu` parameter for the data series in `data`.
   * @param data Data series used to calculate.
   * @return const double The `mu` parameter for the data series.
   */
  const double calculateMu(const std::vector<double>& data) const;
  /**
   * @brief Calculates the `sigma` parameter for the data series in `data`.
   * @param data Data series used to calculate.
   * @param mu The `mu` parameter for the data series.
   * @return const double The `sigma` parameter for the data series.
   */
  const double
  calculateSigma(const std::vector<double>& data, const double& mu) const;
  /**
   * @brief Calculates the maxmimum likelihood values for all model
   * parameters.
   *
   * @param data Data series used to calculate.
   * @return const std::unordered_map<std::string, const double> The maximum
   * likelihood value for model parameters.
   */
  const double
  calculateConditionalVariance(const double& sigma, const double& mu) const;

  /**
   * @brief Calculates the likelihood equation components from data series
   * and model parameters.
   *
   * @param data Data series used to calculate.
   * @param params The model parameters calculated from the data series.
   * @return const GeneralLinearLikelihoodComponents The likelihood equation
   * components.
   */
  const GeneralLinearLikelihoodComponents calculateComponents(
      const std::vector<double>& data, const GeneralLinearParameters& params
  ) const;

  const GeneralLinearParameters
  calculateParameters(const std::vector<double>& data) const;
};

#endif // STOCHASTIC_MODELS_LIKELIHOOD_GENERAL_LINEAR_LIKELIHOOD_H
