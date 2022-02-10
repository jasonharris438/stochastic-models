#ifndef _ORNSTEIN_UHLENBECK_H_
#define _ORNSTEIN_UHLENBECK_H_
#include "stochastic_model.h"

class OrnsteinUhlenbeckModel: public StochasticModel {
private:
    const double mu;
    const double alpha;
    const double sigma;

public:
    virtual const double CoreEquation(const double x, const unsigned int t = 1) override;
};

#endif // _ORNSTEIN_UHLENBECK_H_
