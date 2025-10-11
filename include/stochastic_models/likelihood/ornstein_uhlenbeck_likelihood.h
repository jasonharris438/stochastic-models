#ifndef ORNSTEIN_UHLENBECK_LIKELIHOOD_H
#define ORNSTEIN_UHLENBECK_LIKELIHOOD_H
#include "stochastic_models/likelihood/likelihood.h"
/**
 * @brief Handles calculating maxmimum likelihood parameters specifically under
 * Ornstein-Uhlenbeck processes assumptions for a data series.
 *
 */
class OrnsteinUhlenbeckLikelihood : public LikelihoodCalculatuor {
private:
  /**
   * @brief Calculates likelihood equation components from data parameters and
   * stores in components member.
   *
   * @param data Data series used to calculate maximum likelihood parameters.
   */
  void presetComponents(const std::vector<double>& data);
  /**
   * @brief Calculates likelihood component Xy.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateXy(const std::vector<double>& data) const;
  /**
   * @brief Calculates likelihood component Xx.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateXx(const std::vector<double>& data) const;
  /**
   * @brief Calculates likelihood component Xx.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateXyy(const std::vector<double>& data) const;
  /**
   * @brief Calculates likelihood component Xyy.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateXxx(const std::vector<double>& data) const;
  /**
   * @brief Calculates likelihood component Xxx.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateXxy(const std::vector<double>& data) const;

public:
  OrnsteinUhlenbeckLikelihood();
  /**
   * @brief Calculates the maxmimum likelihood value for the model parameter
   * mu using the values in the components member.
   *
   * @param data Data series used to calculate.
   * @return const double The value for mu parameter.
   */
  const double calculateMu() const;

  /**
   * @brief Calculates the maxmimum likelihood value for the model parameter
   * alpha using the values in the components member and the existing mu
   * value.
   *
   * @param mu The maximum likelihood parmater mu for the model.
   * @return const double The value for alpha parameter.
   */
  const double calculateAlpha(const double& mu) const;
  /**
   * @brief Calculates the maxmimum likelihood value for the model parameter
   * sigma using the values in the components member and the existing mu and
   * alpha values.
   *
   * @param mu The maximum likelihood parmater mu for the model.
   * @param alpha The maximum likelihood parmater alpha for the model.
   * @return const double The value for sigma parameter.
   */
  const double calculateSigma(const double& mu, const double& alpha) const;
  /**
   * @brief Calculates maximum likelihood values for all model parameters.
   *
   * @return const std::unordered_map<std::string, const double> unordered_map
   * containing maximum likelihood model parmeters.
   */
  const std::unordered_map<std::string, const double>
  calculate(const std::vector<double>& data) override;
};
#endif // ORNSTEIN_UHLENBECK_LIKELIHOOD_H
