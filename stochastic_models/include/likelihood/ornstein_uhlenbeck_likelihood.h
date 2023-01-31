#ifndef _ORNSTEIN_UHLENBECK_LIKELIHOOD_H_
#define _ORNSTEIN_UHLENBECK_LIKELIHOOD_H_
#include "likelihood.h"
class OrnsteinUhlenbeckLikelihood : public LikelihoodCalculatuor {
   private:
    const double calculate_x_y(const std::vector<double>& data) const;
    const double calculate_x_x(const std::vector<double>& data) const;
    const double calculate_x_yy(const std::vector<double>& data) const;
    const double calculate_x_xx(const std::vector<double>& data) const;
    const double calculate_x_xy(const std::vector<double>& data) const;
    const double calculate_mu(const std::vector<double>& data) const;

   public:
    OrnsteinUhlenbeckLikelihood();
    const std::map<std::string, double> calculate(
        const std::vector<double>& data) const override;
};
#endif  // _ORNSTEIN_UHLENBECK_LIKELIHOOD_H_
