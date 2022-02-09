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
    const GaussianDistribution *dist;
public:
    StochasticModel();
    ~StochasticModel();
    std::vector<const double> *Simulate(const unsigned int size);
};

class OrnsteinUhlenbeckModel: public StochasticModel {

};

#endif