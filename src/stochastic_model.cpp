#include "stochastic_model.h"
#include "gaussian.h"

#include <vector>

StochasticModel::StochasticModel() {
    dist = new GaussianDistribution();
}

StochasticModel::~StochasticModel(){
    delete dist;
    dist = nullptr;
}

std::vector<double> StochasticModel::Simulate(const double start, const unsigned int size){
    const std::vector<double> distribution_draws = (*dist).Sample(size);
    std::vector<double> vec = {start};

    for(unsigned int n{} ; n < size ; n++){
        const double sample = CoreEquation(vec[n], distribution_draws[n]);
        vec.push_back(sample);
    }

    return vec;
}
