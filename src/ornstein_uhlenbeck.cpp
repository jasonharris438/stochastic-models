#include "stochastic_model.h"
#include "ornstein_uhlenbeck.h"

#include <cmath>

OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel(const double mu, const double alpha, const double sigma) :
    mu(mu), alpha(alpha), sigma(sigma) {}

OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel() :
   OrnsteinUhlenbeckModel(0, 0.01, 0.05) {}

const double OrnsteinUhlenbeckModel::CoreEquation(const double x, const unsigned int t){
    double solution{};
    double delta{std::exp(-alpha * t)};
    std::mt19937 gen;

    const double solution = (x * delta) + (mu * (1 - delta)) + (t * sigma * (*dist).Sample(gen));
    return solution;
}
