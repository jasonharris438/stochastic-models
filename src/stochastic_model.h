#ifndef STOCHASTIC_MODEL_H_
#define STOCHASTIC_MODEL_H_

#include "gaussian.h"

#include <vector>

/**
 * Stochastic Model base class that handles functionality for fitting, analysing,
 * and simulating statistical models.
 */
class StochasticModel {
private:
    GaussianDistribution *dist;
public:
    StochasticModel();
    virtual ~StochasticModel();
    virtual std::vector<float> Simulate(unsigned int size);
};

class OrnsteinUhlenbeckModel: public StochasticModel {

};

#endif