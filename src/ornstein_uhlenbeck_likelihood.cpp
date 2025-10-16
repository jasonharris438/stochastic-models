#include "stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h"

#include "stochastic_models/numeric_utils/helpers.h"

#include <cmath>
#include <iterator>
#include <new>
#include <numeric>

OrnsteinUhlenbeckLikelihoodComponents::OrnsteinUhlenbeckLikelihoodComponents(
    const double lead_sum,
    const double lag_sum,
    const double lead_sum_squared,
    const double lag_sum_squared,
    const double lead_lag_sum_product,
    const uint32_t n_obs
)
    : lead_sum(lead_sum), lag_sum(lag_sum), lead_sum_squared(lead_sum_squared),
      lag_sum_squared(lag_sum_squared),
      lead_lag_sum_product(lead_lag_sum_product), n_obs(n_obs) {}

const double OrnsteinUhlenbeckLikelihood::calculateLeadSum(
    const std::vector<double>& data
) const {
  std::vector<double>::const_iterator iter = data.cbegin();
  std::advance(iter, 1);
  return std::reduce(iter, data.cend(), 0.0);
};
const double OrnsteinUhlenbeckLikelihood::calculateLagSum(
    const std::vector<double>& data
) const {
  std::vector<double>::const_iterator iter = data.cend();
  std::advance(iter, -1);
  return std::reduce(data.cbegin(), iter, 0.0);
};
const double OrnsteinUhlenbeckLikelihood::calculateLeadSumSquared(
    const std::vector<double>& data
) const {
  const std::vector<double> squared = valuesSquared(data);
  std::vector<double>::const_iterator iter = squared.cbegin();
  std::advance(iter, 1);
  return std::reduce(iter, squared.cend(), 0.0);
};
const double OrnsteinUhlenbeckLikelihood::calculateLagSumSquared(
    const std::vector<double>& data
) const {
  const std::vector<double> squared = valuesSquared(data);
  std::vector<double>::const_iterator iter = squared.cend();
  std::advance(iter, -1);
  return std::reduce(squared.cbegin(), iter, 0.0);
};
const double OrnsteinUhlenbeckLikelihood::calculateLeadLagSumProduct(
    const std::vector<double>& data
) const {
  std::vector<double>::const_iterator iter_lead = data.cbegin();
  std::advance(iter_lead, 1);

  std::vector<double>::const_iterator iter_lag_start = data.cbegin();
  std::vector<double>::const_iterator iter_lag_end = data.cend();
  std::advance(iter_lag_end, -1);

  return std::inner_product(iter_lag_start, iter_lag_end, iter_lead, 0.0);
};
OrnsteinUhlenbeckLikelihoodComponents
OrnsteinUhlenbeckLikelihood::calculateComponents(
    const std::vector<double>& data
) {
  return OrnsteinUhlenbeckLikelihoodComponents{
      calculateLeadSum(data),           calculateLagSum(data),
      calculateLeadSumSquared(data),    calculateLagSumSquared(data),
      calculateLeadLagSumProduct(data), static_cast<int>(data.size())
  };
};
const double OrnsteinUhlenbeckLikelihood::calculateMu(
    const OrnsteinUhlenbeckLikelihoodComponents& components
) const {
  return ((components.lead_sum * components.lag_sum_squared) -
          (components.lag_sum * components.lead_lag_sum_product)) /
         ((components.n_obs *
           (components.lag_sum_squared - components.lead_lag_sum_product)) -
          (std::pow(components.lag_sum, 2) -
           (components.lead_sum * components.lag_sum)));
};
const double OrnsteinUhlenbeckLikelihood::calculateAlpha(
    const OrnsteinUhlenbeckLikelihoodComponents& components, const double& mu
) const {
  return std::log(
             components.lag_sum_squared - (2 * mu * components.lag_sum) +
             (components.n_obs * (std::pow(mu, 2)))
         ) -
         std::log(
             components.lead_lag_sum_product - (mu * components.lag_sum) -
             (mu * components.lead_sum) + (components.n_obs * (std::pow(mu, 2)))
         );
};
const double OrnsteinUhlenbeckLikelihood::calculateSigma(
    const OrnsteinUhlenbeckLikelihoodComponents& components,
    const double& mu,
    const double& alpha
) const {
  const double exp_alpha{std::exp(-alpha)};
  double sigma{
      components.lead_sum_squared -
      (2 * exp_alpha * components.lead_lag_sum_product) +
      ((std::pow(exp_alpha, 2)) * components.lag_sum_squared) -
      ((2 * mu * (1 - exp_alpha)) *
       (components.lead_sum - (exp_alpha * components.lag_sum))) +
      (components.n_obs * std::pow(mu, 2) * std::pow(1 - exp_alpha, 2))
  };
  sigma *= (1 / components.n_obs);
  sigma *= ((2 * exp_alpha) / (1 - std::pow(exp_alpha, 2)));
  return sigma;
};
const OrnsteinUhlenbeckParameters
calculateOrnsteinUhlenbeckParameters(const std::vector<double>& data) {
  OrnsteinUhlenbeckLikelihood ou_likelihood;
  OrnsteinUhlenbeckLikelihoodComponents components =
      ou_likelihood.calculateComponents(data);
  const double mu = ou_likelihood.calculateMu(components);
  const double alpha = ou_likelihood.calculateAlpha(components, mu);
  const double sigma = ou_likelihood.calculateSigma(components, mu, alpha);
  return OrnsteinUhlenbeckParameters{mu, alpha, sigma};
};
