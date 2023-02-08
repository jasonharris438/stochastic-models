#ifndef _ORNSTEIN_UHLENBECK_LIKELIHOOD_H_
#define _ORNSTEIN_UHLENBECK_LIKELIHOOD_H_
#include "likelihood.h"
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
    void preset_components(const std::vector<double>& data);
    /**
     * @brief Calculates likelihood component x_y.
     *
     * @param data Data series used to calculate.
     * @return const double Likelihood component value.
     */
    const double calculate_x_y(const std::vector<double>& data) const;
    /**
     * @brief Calculates likelihood component x_x.
     *
     * @param data Data series used to calculate.
     * @return const double Likelihood component value.
     */
    const double calculate_x_x(const std::vector<double>& data) const;
    /**
     * @brief Calculates likelihood component x_x.
     *
     * @param data Data series used to calculate.
     * @return const double Likelihood component value.
     */
    const double calculate_x_yy(const std::vector<double>& data) const;
    /**
     * @brief Calculates likelihood component x_yy.
     *
     * @param data Data series used to calculate.
     * @return const double Likelihood component value.
     */
    const double calculate_x_xx(const std::vector<double>& data) const;
    /**
     * @brief Calculates likelihood component x_xx.
     *
     * @param data Data series used to calculate.
     * @return const double Likelihood component value.
     */
    const double calculate_x_xy(const std::vector<double>& data) const;

   public:
    OrnsteinUhlenbeckLikelihood();
    /**
     * @brief Calculates the maxmimum likelihood value for the model parameter
     * mu using the values in the components member.
     *
     * @param data Data series used to calculate.
     * @return const double The value for mu parameter.
     */
    const double calculate_mu() const;

    /**
     * @brief Calculates the maxmimum likelihood value for the model parameter
     * alpha using the values in the components member and the existing mu
     * value.
     *
     * @param mu The maximum likelihood parmater mu for the model.
     * @return const double The value for alpha parameter.
     */
    const double calculate_alpha(const double& mu) const;
    /**
     * @brief Calculates the maxmimum likelihood value for the model parameter
     * sigma using the values in the components member and the existing mu and
     * alpha values.
     *
     * @param mu The maximum likelihood parmater mu for the model.
     * @param alpha The maximum likelihood parmater alpha for the model.
     * @return const double The value for sigma parameter.
     */
    const double calculate_sigma(const double& mu, const double& alpha) const;
    /**
     * @brief Calculates maximum likelihood values for all model parameters.
     *
     * @return const std::unordered_map<std::string, const double> unordered_map
     * containing maximum likelihood model parmeters.
     */
    const std::unordered_map<std::string, const double> calculate(
        const std::vector<double>& data) override;
};
#endif  // _ORNSTEIN_UHLENBECK_LIKELIHOOD_H_
