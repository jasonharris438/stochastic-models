#ifndef STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_LIKELIHOOD_H
#define STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_LIKELIHOOD_H

#include <cstdint>
#include <vector>

struct OrnsteinUhlenbeckParameters {
  double mu;
  double alpha;
  double sigma;
};

struct OrnsteinUhlenbeckLikelihoodComponents {
  double lead_sum;
  double lag_sum;
  double lead_sum_squared;
  double lag_sum_squared;
  double lead_lag_sum_product;
  uint32_t n_obs;
};

class OrnsteinUhlenbeckLikelihoodComponentCalculator {
public:
  /**
   * @brief Calculates the sum of the lead series used in likelihood
   * calculations.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateLeadSum(const std::vector<double>& data) const;
  /**
   * @brief Calculates the sum of the lagged series used in likelihood
   * calculations.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateLagSum(const std::vector<double>& data) const;
  /**
   * @brief Calculates the sum of squares of the lead series used in likelihood
   * calculations.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateLeadSumSquared(const std::vector<double>& data) const;
  /**
   * @brief Calculates the sum of squares of the lagged series used in
   * likelihood calculations.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateLagSumSquared(const std::vector<double>& data) const;
  /**
   * @brief Calculates the sum of the products of the lead and lagged series
   * used in likelihood calculations.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double
  calculateLeadLagSumProduct(const std::vector<double>& data) const;
  const double
  updateLeadSum(const double& lead_sum, const double& new_observation) const;
  const double
  updateLagSum(const double& lag_sum, const double& last_observation) const;
  const double updateLeadSumSquared(
      const double& lead_sum_squared, const double& new_observation
  ) const;
  const double updateLagSumSquared(
      const double& lag_sum_squared, const double& last_observation
  ) const;
  const double updateLeadLagSumProduct(
      const double& lead_lag_sum_product,
      const double& new_observation,
      const double& last_observation
  ) const;
  /**
   * @brief Calculates the maxmimum likelihood value for the model parameter
   * mu using the values in the components member.
   *
   * @param components Summary values of a data series used to calculate the
   * maximum likelihood parameters.
   * @return const double The value for mu parameter.
   */
  const double
  calculateMu(const OrnsteinUhlenbeckLikelihoodComponents& components) const;

  /**
   * @brief Calculates the maxmimum likelihood value for the model parameter
   * alpha using the values in the components member and the existing mu
   * value.
   *
   * @param components Summary values of a data series used to calculate the
   * maximum likelihood parameters.
   * @param mu The maximum likelihood parmater mu for the model.
   * @return const double The value for alpha parameter.
   */
  const double calculateAlpha(
      const OrnsteinUhlenbeckLikelihoodComponents& components, const double& mu
  ) const;
  /**
   * @brief Calculates the maxmimum likelihood value for the model parameter
   * sigma using the values in the components member and the existing mu and
   * alpha values.
   *
   * @param components Summary values of a data series used to calculate the
   * maximum likelihood parameters.
   * @param mu The maximum likelihood parmater mu for the model.
   * @param alpha The maximum likelihood parmater alpha for the model.
   * @return const double The value for sigma parameter.
   */
  const double calculateSigma(
      const OrnsteinUhlenbeckLikelihoodComponents& components,
      const double& mu,
      const double& alpha
  ) const;
};

/**
 * @brief Handles calculating maxmimum likelihood parameters specifically under
 * Ornstein-Uhlenbeck processes assumptions for a data series.
 *
 */
class OrnsteinUhlenbeckLikelihood {
private:
  /**
   * @brief Component calculator instance to use when calculating likelihood
   * components from a data series.
   */
  OrnsteinUhlenbeckLikelihoodComponentCalculator component_calculator;

public:
  /**
   * @brief Calculates the internal likelihood equation components from a data
   * series.
   *
   * @param data Data series used to calculate maximum likelihood parameters.
   */
  OrnsteinUhlenbeckLikelihoodComponents
  calculateComponents(const std::vector<double>& data) const;

  const OrnsteinUhlenbeckParameters calculateParameters(
      const OrnsteinUhlenbeckLikelihoodComponents& components
  ) const;
  const OrnsteinUhlenbeckLikelihoodComponents updateComponents(
      const OrnsteinUhlenbeckLikelihoodComponents& components,
      const double& new_observation,
      const double& last_observation
  ) const;
};

#endif // STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_LIKELIHOOD_H
