#include "likelihood/ornstein_uhlenbeck_likelihood.h"

#include <cmath>
#include <iterator>
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
const double OrnsteinUhlenbeckLikelihood::calculate_mu(
    const std::vector<double>& data) const {
    const double x_y{calculate_x_y(data)};
    const double x_x{calculate_x_x(data)};
    const double x_xx{calculate_x_xx(data)};
    const double x_yy{calculate_x_yy(data)};
    const double x_xy{calculate_x_xy(data)};
    const int size_of{data.size()};

    return ((x_y * x_xx) - (x_x * x_xy)) /
           ((size_of * (x_xx - x_xy)) - (std::pow(x_x, 2) - (x_y * x_x)));
};
const std::map<std::string, double> OrnsteinUhlenbeckLikelihood::calculate(
    const std::vector<double>& data) const {
    const double mu{calculate_mu(data)};
    std::map<std::string, double> m{{"mu", mu}};
    return m;
};
