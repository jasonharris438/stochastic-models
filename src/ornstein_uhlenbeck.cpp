#include "stochastic_models/sde/ornstein_uhlenbeck.h"

#include <cmath>

#include "stochastic_models/distributions/gaussian.h"
/**
 * @brief No args constructor delegates to main constructor.
 *
 */
OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel()
    : OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel(0.0, 0.0, 1.0) {}
/**
 * @brief Main constructor.
 *
 */
OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel(const double mu,
                                               const double alpha,
                                               const double sigma)
    : mu(mu), alpha(alpha), sigma(sigma) {
    dist = new GaussianDistribution();
}
/**
 * @brief Copy constructor.
 *
 * @param other The OrnsteinUhlenbeckModel to copy.
 */
OrnsteinUhlenbeckModel::OrnsteinUhlenbeckModel(
    const OrnsteinUhlenbeckModel& other)
    : mu(other.mu), alpha(other.alpha), sigma(other.sigma) {
    dist = new GaussianDistribution();
}
/**
 * @brief Destructor.
 *
 */
OrnsteinUhlenbeckModel::~OrnsteinUhlenbeckModel() {
    delete dist;
    dist = nullptr;
}
const OrnsteinUhlenbeckModel* OrnsteinUhlenbeckModel::clone() const {
    return new OrnsteinUhlenbeckModel(*this);
}
const double OrnsteinUhlenbeckModel::getMean() const { return mu; }
const double OrnsteinUhlenbeckModel::getUnconditionalVariance() const {
    return std::pow(sigma, 2) / (2 * alpha);
}
std::vector<double> OrnsteinUhlenbeckModel::Simulate(
    const double start, const unsigned int& size, const unsigned int& t) const {
    const std::vector<double> distribution_draws = (*dist).sample(size - 1);
    std::vector<double> vec = {start};

    for (const unsigned int& val : distribution_draws) {
        const double& last = vec.back();
        const double sample = coreEquation(last, val, t);
        vec.push_back(sample);
    }

    return vec;
}
const double OrnsteinUhlenbeckModel::coreEquation(const double& x,
                                                  const double& noise,
                                                  const unsigned int& t) const {
    double delta{std::exp(-alpha * t)};
    const double solution =
        (x * delta) + (mu * (1 - delta)) + (t * sigma * noise);
    return solution;
}
