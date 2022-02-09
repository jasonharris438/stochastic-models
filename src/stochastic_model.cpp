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

std::vector<const double> *StochasticModel::Simulate(const unsigned int size){
    // Dynamic heap allocation for vector
    std::vector<const double> *vec_ptr {nullptr};
    vec_ptr = new std::vector<const double>;

    // Seed random number generator with current time
    long int number = static_cast<long int> (time(NULL));
    std::mt19937 gen(number);

    for(unsigned int n{} ; n < size ; n++){
        const double sample = (*dist).Sample(gen);
        (*vec_ptr).push_back(sample);
    }

    return vec_ptr;
}