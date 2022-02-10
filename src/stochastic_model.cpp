#include "stochastic_model.h"
#include "gaussian.h"

#include <vector>
#include <time.h>

StochasticModel::StochasticModel() {
    dist = new GaussianDistribution();
}

StochasticModel::~StochasticModel(){
    delete dist;
    dist = nullptr;
}

const double StochasticModel::Simulate(){
    // Seed random number generator with current time
    long int number = static_cast<long int> (time(NULL));
    std::mt19937 gen(number);

    const double sample = (*dist).Sample(gen);
    return sample;
}
