#include "stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h"

#include <cmath>
#include <iterator>
#include <new>
#include <numeric>

#include "stochastic_models/numeric_utils/helpers.h"
OrnsteinUhlenbeckLikelihood::OrnsteinUhlenbeckLikelihood(){};
const double OrnsteinUhlenbeckLikelihood::calculateXy(
    const std::vector<double>& data) const {
    std::vector<double>::const_iterator iter = data.cbegin();
    std::advance(iter, 1);
    return std::reduce(iter, data.cend(), 0.0);
};
const double OrnsteinUhlenbeckLikelihood::calculateXx(
    const std::vector<double>& data) const {
    std::vector<double>::const_iterator iter = data.cend();
    std::advance(iter, -1);
    return std::reduce(data.cbegin(), iter, 0.0);
};
const double OrnsteinUhlenbeckLikelihood::calculateXyy(
    const std::vector<double>& data) const {
    std::vector<double> new_vec(data);
    valuesSquared(new_vec);
    std::vector<double>::const_iterator iter = new_vec.cbegin();
    std::advance(iter, 1);
    return std::reduce(iter, new_vec.cend(), 0.0);
};
const double OrnsteinUhlenbeckLikelihood::calculateXxx(
    const std::vector<double>& data) const {
    std::vector<double> new_vec(data);
    valuesSquared(new_vec);
    std::vector<double>::const_iterator iter = new_vec.cend();
    std::advance(iter, -1);
    return std::reduce(new_vec.cbegin(), iter, 0.0);
};
const double OrnsteinUhlenbeckLikelihood::calculateXxy(
    const std::vector<double>& data) const {
    std::vector<double>::const_iterator iter_y = data.cbegin();
    std::advance(iter_y, 1);

    std::vector<double>::const_iterator iter_x_start = data.cbegin();
    std::vector<double>::const_iterator iter_x_end = data.cend();
    std::advance(iter_x_end, -1);

    return std::inner_product(iter_x_start, iter_x_end, iter_y, 0.0);
};
void OrnsteinUhlenbeckLikelihood::presetComponents(
    const std::vector<double>& data) {
    if (components == nullptr) {
        components = new const std::unordered_map<std::string, const double>{
            {"x_y", calculateXy(data)},   {"x_x", calculateXx(data)},
            {"x_yy", calculateXyy(data)}, {"x_xx", calculateXxx(data)},
            {"x_xy", calculateXxy(data)}, {"n", data.size()}};
    } else {
        throw std::bad_alloc();
    }
};
const double OrnsteinUhlenbeckLikelihood::calculateMu() const {
    return ((components->at("x_y") * components->at("x_xx")) -
            (components->at("x_x") * components->at("x_xy"))) /
           ((components->at("n") *
             (components->at("x_xx") - components->at("x_xy"))) -
            (std::pow(components->at("x_x"), 2) -
             (components->at("x_y") * components->at("x_x"))));
};
const double OrnsteinUhlenbeckLikelihood::calculateAlpha(
    const double& mu) const {
    return std::log(components->at("x_xx") - (2 * mu * components->at("x_x")) +
                    (components->at("n") * (std::pow(mu, 2)))) -
           std::log(components->at("x_xy") - (mu * components->at("x_x")) -
                    (mu * components->at("x_y")) +
                    (components->at("n") * (std::pow(mu, 2))));
};
const double OrnsteinUhlenbeckLikelihood::calculateSigma(
    const double& mu, const double& alpha) const {
    const double exp_alpha{std::exp(-alpha)};
    double sigma{
        components->at("x_yy") - (2 * exp_alpha * components->at("x_xy")) +
        ((std::pow(exp_alpha, 2)) * components->at("x_xx")) -
        ((2 * mu * (1 - exp_alpha)) *
         (components->at("x_y") - (exp_alpha * components->at("x_x")))) +
        (components->at("n") * std::pow(mu, 2) * std::pow(1 - exp_alpha, 2))};
    sigma *= (1 / components->at("n"));
    sigma *= ((2 * exp_alpha) / (1 - std::pow(exp_alpha, 2)));
    return sigma;
};
const std::unordered_map<std::string, const double>
OrnsteinUhlenbeckLikelihood::calculate(const std::vector<double>& data) {
    presetComponents(data);
    const double mu{calculateMu()};
    const double alpha{calculateAlpha(mu)};
    const double sigma{calculateSigma(mu, alpha)};
    std::unordered_map<std::string, const double> m{
        {"mu", mu}, {"alpha", alpha}, {"sigma", sigma}};
    return m;
};
