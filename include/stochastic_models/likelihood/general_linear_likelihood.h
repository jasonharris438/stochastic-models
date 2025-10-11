#ifndef GENERAL_LINEAR_LIKELIHOOD_H
#define GENERAL_LINEAR_LIKELIHOOD_H
#include "stochastic_models/likelihood/likelihood.h"

/**
 * @brief Handles calculating maxmimum likelihood parameters specifically under
 * General Linear SDE process assumptions for a data series.
 *
 */
class GeneralLinearLikelihood : public LikelihoodCalculatuor {
public:
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
   * @brief Calculates the kernel value for the sigma parameter of the data
   * series in `data`.
   * @param n_observations The number of observations in the data series.
   * @param sigma The sigma parameter for the data series.
   * @return const double The kernel value for the sigma parameter.
   */
  const double calculateSigmaKernel(
      const u_int32_t& n_observations, const double& sigma
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
  const std::unordered_map<std::string, const double>
  calculate(const std::vector<double>& data) override;
  const double
  calculateConditionalVariance(const double& sigma, const double& mu) const;
};
#endif // GENERAL_LINEAR_LIKELIHOOD_H
