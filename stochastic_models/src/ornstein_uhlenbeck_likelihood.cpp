#include "likelihood/ornstein_uhlenbeck_likelihood.h"

#include <iterator>
#include <numeric>

OrnsteinUhlenbeckLikelihood::OrnsteinUhlenbeckLikelihood(){};
double OrnsteinUhlenbeckLikelihood::calculate_x_y(
    const std::vector<double>& data) const {
    std::vector<double>::const_iterator iter = data.cbegin();
    advance(iter, 1);
    return std::reduce(iter, data.cend(), 0.0);
};
double OrnsteinUhlenbeckLikelihood::calculate_x_x(
    const std::vector<double>& data) const {
    std::vector<double>::const_iterator iter = data.cend();
    advance(iter, -1);
    return std::reduce(data.cbegin(), iter, 0.0);
};
const std::map<std::string, double> OrnsteinUhlenbeckLikelihood::calculate(
    const std::vector<double>& data) const {
    const double x_y{OrnsteinUhlenbeckLikelihood::calculate_x_y(data)};
    const double x_x{OrnsteinUhlenbeckLikelihood::calculate_x_x(data)};
    std::map<std::string, double> m{{"x", x_y + x_x}, {"GPU", 15.1}};
    return m;
};
