#include "likelihood/ornstein_uhlenbeck_likelihood.h"

#include <cmath>
#include <iterator>
#include <new>
#include <numeric>

#include "numeric_utils/helpers.h"
OrnsteinUhlenbeckLikelihood::OrnsteinUhlenbeckLikelihood(){};
const double OrnsteinUhlenbeckLikelihood::calculate_x_y(
    const std::vector<double>& data) const {
    std::vector<double>::const_iterator iter = data.cbegin();
    advance(iter, 1);
    return std::reduce(iter, data.cend(), 0.0);
};
const double OrnsteinUhlenbeckLikelihood::calculate_x_x(
    const std::vector<double>& data) const {
    std::vector<double>::const_iterator iter = data.cend();
    advance(iter, -1);
    return std::reduce(data.cbegin(), iter, 0.0);
};
const double OrnsteinUhlenbeckLikelihood::calculate_x_yy(
    const std::vector<double>& data) const {
    std::vector<double> new_vec(data);
    ValuesSquared(new_vec);
    std::vector<double>::const_iterator iter = new_vec.cbegin();
    advance(iter, 1);
    return std::reduce(iter, new_vec.cend(), 0.0);
};
const double OrnsteinUhlenbeckLikelihood::calculate_x_xx(
    const std::vector<double>& data) const {
    std::vector<double> new_vec(data);
    ValuesSquared(new_vec);
    std::vector<double>::const_iterator iter = new_vec.cend();
    advance(iter, -1);
    return std::reduce(new_vec.cbegin(), iter, 0.0);
};
const double OrnsteinUhlenbeckLikelihood::calculate_x_xy(
    const std::vector<double>& data) const {
    std::vector<double> new_vec(data);

    std::vector<double>::const_iterator iter_y = new_vec.cbegin();
    advance(iter_y, 1);

    std::vector<double>::const_iterator iter_x_start = new_vec.cbegin();
    std::vector<double>::const_iterator iter_x_end = new_vec.cend();
    advance(iter_x_end, -1);

    return std::inner_product(iter_x_start, iter_x_end, iter_y, 0.0);
};
void OrnsteinUhlenbeckLikelihood::preset_components(
    const std::vector<double>& data) {
    if (components == nullptr) {
        components = new const std::map<std::string, const double>{
            {"x_y", calculate_x_y(data)},   {"x_x", calculate_x_x(data)},
            {"x_yy", calculate_x_yy(data)}, {"x_xx", calculate_x_xx(data)},
            {"x_xy", calculate_x_xy(data)}, {"n", data.size()}};
    } else {
        throw std::bad_alloc();
    }
};
const double OrnsteinUhlenbeckLikelihood::calculate_mu() const {
    return ((components->at("x_y") * components->at("x_xx")) -
            (components->at("x_x") * components->at("x_xy"))) /
           ((components->at("n") *
             (components->at("x_xx") - components->at("x_xy"))) -
            (std::pow(components->at("x_x"), 2) -
             (components->at("x_y") * components->at("x_x"))));
};
const double OrnsteinUhlenbeckLikelihood::calculate_alpha(
    const double& mu) const {
    return std::log(components->at("x_xx") - (2 * mu * components->at("x_x")) +
                    (components->at("n") * (std::pow(mu, 2)))) -
           std::log(components->at("x_xy") - (mu * components->at("x_x")) -
                    (mu * components->at("x_y")) +
                    (components->at("n") * (std::pow(mu, 2))));
};
const double OrnsteinUhlenbeckLikelihood::calculate_sigma(
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
const std::map<std::string, const double>
OrnsteinUhlenbeckLikelihood::calculate(const std::vector<double>& data) {
    preset_components(data);
    const double mu{calculate_mu()};
    const double alpha{calculate_alpha(mu)};
    const double sigma{calculate_sigma(mu, alpha)};
    std::map<std::string, const double> m{
        {"mu", mu}, {"alpha", alpha}, {"sigma", sigma}};
    return m;
};
