#include "likelihood/ornstein_uhlenbeck_likelihood.h"

#include <numeric>

OrnsteinUhlenbeckLikelihood::OrnsteinUhlenbeckLikelihood(){};
double OrnsteinUhlenbeckLikelihood::calculate_x_y(
    const std::vector<double>& data) const {
    return std::reduce(data.begin() + 1, data.end(), 0, std::plus<double>());
};
const std::map<std::string, double>& OrnsteinUhlenbeckLikelihood::calculate(
    const std::vector<double>& data) const {
    std::map<std::string, double> m{{"CPU", 10.0}, {"GPU", 15.1}};
    return m;
    ;
};
