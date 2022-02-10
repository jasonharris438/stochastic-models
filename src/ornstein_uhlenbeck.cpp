#include "stochastic_model.h"
#include "ornstein_uhlenbeck.h"

#include <cmath>

double OrnsteinUhlenbeckModel::CoreEquation(const double x, const unsigned int t){
    double solution{};
    double delta{std::exp(-alpha * t)};
    std::mt19937 gen;

    solution = (x * delta) + (mu * (1 - delta)) + (sigma * (*dist).Sample(gen));
    return solution;
}
