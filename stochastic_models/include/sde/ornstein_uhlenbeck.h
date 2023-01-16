#ifndef _ORNSTEIN_UHLENBECK_H_
#define _ORNSTEIN_UHLENBECK_H_
#include "sde/stochastic_model.h"

class OrnsteinUhlenbeckModel : public StochasticModel {
   private:
    const double mu;
    const double alpha;
    const double sigma;

   public:
    OrnsteinUhlenbeckModel();
    OrnsteinUhlenbeckModel(const double mu, const double alpha,
                           const double sigma);
    const double CoreEquation(const double x, const double noise,
                              const unsigned int t = 1) const override;
    std::vector<double> Simulate(const double start,
                                 const unsigned int size) const override;
    ~OrnsteinUhlenbeckModel();
};

#endif  // _ORNSTEIN_UHLENBECK_H_
