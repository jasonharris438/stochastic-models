#include "stochastic_model.h"
#include "ornstein_uhlenbeck.h"

#include <cmath>

const double OrnsteinUhlenbeckModel::CoreEquation(const double x, const unsigned int t){
    double solution{};
    double delta{std::exp(-alpha * t)};

    const double solution = (x * delta) + (mu * (1 - delta)) + (t * sigma * Simulate());
    return solution;
}
