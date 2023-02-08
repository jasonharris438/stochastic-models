#include "sde/stochastic_model.h"

StochasticModel::~StochasticModel() {
    delete dist;
    dist = nullptr;
};
