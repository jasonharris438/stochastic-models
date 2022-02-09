#include <vector>
#include "stochastic_model.h"
#include "gaussian.h"


StochasticModel::StochasticModel() {
    dist = new GaussianDistribution();
}

StochasticModel::~StochasticModel(){
        delete dist;
        dist = nullptr;
}
