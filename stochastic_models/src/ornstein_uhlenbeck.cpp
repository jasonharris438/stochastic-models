#include "sde/ornstein_uhlenbeck.h"

#include <cmath>

#include "distributions/gaussian.h"
OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel()
    : OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel(0, 0, 1) {}
OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel(const double mu,
                                               const double alpha,
                                               const double sigma)
    : mu(mu), alpha(alpha), sigma(sigma) {
    dist = new GaussianDistribution();
}
OrnsteinUhlenbeckModel::~OrnsteinUhlenbeckModel() {}
std::vector<double> OrnsteinUhlenbeckModel::Simulate(
    const double start, const unsigned int size) const {
    const std::vector<double> distribution_draws = (*dist).Sample(size);
    std::vector<double> vec = {start};

    for (unsigned int n{}; n < size; n++) {
        const double sample = CoreEquation(vec[n], distribution_draws[n]);
        vec.push_back(sample);
    }

    return vec;
}
const double OrnsteinUhlenbeckModel::CoreEquation(const double x,
                                                  const double noise,
                                                  const unsigned int t) const {
    double delta{std::exp(-alpha * t)};
    const double solution =
        (x * delta) + (mu * (1 - delta)) + (t * sigma * noise);
    return solution;
}
