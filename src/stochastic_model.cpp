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

std::vector<double> StochasticModel::Simulate(const unsigned int size){
    const double distribution_mean = (*dist).getMean();
    std::vector<double> vec = {distribution_mean};

    for(unsigned int n{} ; n < size ; n++){
        const double sample = CoreEquation(vec[n]);
        vec.push_back(sample);
    }

    return vec;
}
