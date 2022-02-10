#ifndef _STOCHASTIC_MODEL_H_
#define _STOCHASTIC_MODEL_H_
#include "gaussian.h"

#include <vector>

/**
 * Stochastic Model base class that handles functionality for fitting, analysing,
 * and simulating statistical models.
 */
class StochasticModel {
private:
    const GaussianDistribution *dist;
public:
    StochasticModel();
    std::vector<const double> Simulate(const unsigned int size);
    virtual double CoreEquation(const double x) = 0;
    virtual ~StochasticModel();
};

#endif // _STOCHASTIC_MODEL_H_
