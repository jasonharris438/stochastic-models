#ifndef ORNSTEIN_UHLENBECK_LIKELIHOOD_H
#define ORNSTEIN_UHLENBECK_LIKELIHOOD_H

#include <cstdint>
#include <vector>

struct OrnsteinUhlenbeckParameters {
  double mu;
  double alpha;
  double sigma;
};

class OrnsteinUhlenbeckLikelihoodComponents {
public:
  const double lead_sum;
  const double lag_sum;
  const double lead_sum_squared;
  const double lag_sum_squared;
  const double lead_lag_sum_product;
  const uint32_t n_obs;

  OrnsteinUhlenbeckLikelihoodComponents(
      const double lead_sum,
      const double lag_sum,
      const double lead_sum_squared,
      const double lag_sum_squared,
      const double lead_lag_sum_product,
      const uint32_t n_obs
  );
};

/**
 * @brief Handles calculating maxmimum likelihood parameters specifically under
 * Ornstein-Uhlenbeck processes assumptions for a data series.
 *
 */
class OrnsteinUhlenbeckLikelihood {
private:
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

public:
  /**
   * @brief Calculates likelihood equation components from data parameters and
   * stores in components member.
   *
   * @param data Data series used to calculate maximum likelihood parameters.
   */
  OrnsteinUhlenbeckLikelihoodComponents
  calculateComponents(const std::vector<double>& data);
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

  const OrnsteinUhlenbeckParameters
  calculateParameters(const std::vector<double>& data);
};

#endif // ORNSTEIN_UHLENBECK_LIKELIHOOD_H
