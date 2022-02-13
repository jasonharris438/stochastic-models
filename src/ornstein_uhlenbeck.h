#ifndef _ORNSTEIN_UHLENBECK_H_
#define _ORNSTEIN_UHLENBECK_H_
#include "stochastic_model.h"

class OrnsteinUhlenbeckModel: public StochasticModel {
private:
    const double mu;
    const double alpha;
    const double sigma;

public:
    OrnsteinUhlenbeckModel();
    OrnsteinUhlenbeckModel(const double mu, const double alpha, const double sigma);
    virtual const double CoreEquation(const double x, const double noise, const unsigned int t = 1) override;
    ~OrnsteinUhlenbeckModel();
};

#endif // _ORNSTEIN_UHLENBECK_H_
