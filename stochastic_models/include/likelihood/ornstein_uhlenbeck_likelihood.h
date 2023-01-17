#ifndef _ORNSTEIN_UHLENBECK_LIKELIHOOD_H_
#define _ORNSTEIN_UHLENBECK_LIKELIHOOD_H_
#include "likelihood.h"

class OrnsteinUhlenbeckLikelihood : public LikelihoodCalculatuor {
   private:
    double calculate_x_y(const std::vector<double>& data) const;
    double calculate_x_x(const std::vector<double>& data) const;
    // double calculate_mu()

   public:
    OrnsteinUhlenbeckLikelihood();
    const std::map<std::string, double> calculate(
        const std::vector<double>& data) const override;
};

#endif  // _ORNSTEIN_UHLENBECK_LIKELIHOOD_H_
