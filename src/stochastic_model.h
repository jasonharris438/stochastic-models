#ifndef _STOCHASTIC_MODEL_H_
#define _STOCHASTIC_MODEL_H_
#include "gaussian.h"

/**
 * Stochastic Model base class that handles functionality for fitting, analysing,
 * and simulating statistical models.
 */
class StochasticModel {
private:
    const GaussianDistribution *dist;
public:
    StochasticModel();
    const double Simulate();
    virtual const double CoreEquation(const double x, const unsigned int t = 1) = 0;
    virtual ~StochasticModel();
};

#endif // _STOCHASTIC_MODEL_H_
