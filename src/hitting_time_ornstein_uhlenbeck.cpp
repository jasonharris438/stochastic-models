#include "stochastic_models/hitting_times/hitting_time_ornstein_uhlenbeck.h"

#include <cmath>

HittingTimeOrnsteinUhlenbeck::HittingTimeOrnsteinUhlenbeck(const double mu,
                                                           const double alpha,
                                                           const double sigma)
    : mu(mu), alpha(alpha), sigma(sigma) {}
HittingTimeOrnsteinUhlenbeck::HittingTimeOrnsteinUhlenbeck(
    const HittingTimeOrnsteinUhlenbeck& other)
    : mu(other.mu), alpha(other.alpha), sigma(other.sigma) {}
const HittingTimeOrnsteinUhlenbeck* HittingTimeOrnsteinUhlenbeck::clone()
    const {
    return new HittingTimeOrnsteinUhlenbeck(*this);
}
const double HittingTimeOrnsteinUhlenbeck::hittingTimeDensityCore(
    const double& x) const {
    return exp(x * alpha * (x - 2 * mu) / (pow(sigma, 2)));
}
const double HittingTimeOrnsteinUhlenbeck::optimalTradingFCore(
    const double& x, const double& u, const double& r) const {
    return pow(u, (r / alpha) - 1) *
           exp(sqrt(2 * alpha / pow(sigma, 2)) * (x - mu) * u -
               (pow(u, 2) / 2));
}
const double HittingTimeOrnsteinUhlenbeck::optimalTradingGCore(
    const double& x, const double& u, const double& r) const {
    return pow(u, (r / alpha) - 1) *
           exp(sqrt(2 * alpha / pow(sigma, 2)) * (mu - x) * u -
               (pow(u, 2) / 2));
}
const double HittingTimeOrnsteinUhlenbeck::optimalTradingLCore(
    const double& r, const double& c) const {
    return ((alpha * mu) + (r * c)) / (r + alpha);
}
