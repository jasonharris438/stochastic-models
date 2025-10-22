#include "stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h"

#include "stochastic_models/numeric_utils/helpers.h"

#include <cmath>
#include <iterator>
#include <new>
#include <numeric>

const double OrnsteinUhlenbeckLikelihoodComponentCalculator::calculateLeadSum(
    const std::vector<double>& data
) const {
  std::vector<double>::const_iterator iter = data.cbegin();
  std::advance(iter, 1);
  return std::reduce(iter, data.cend(), 0.0);
};
const double OrnsteinUhlenbeckLikelihoodComponentCalculator::calculateLagSum(
    const std::vector<double>& data
) const {
  std::vector<double>::const_iterator iter = data.cend();
  std::advance(iter, -1);
  return std::reduce(data.cbegin(), iter, 0.0);
};
const double
OrnsteinUhlenbeckLikelihoodComponentCalculator::calculateLeadSumSquared(
    const std::vector<double>& data
) const {
  const std::vector<double> squared = valuesSquared(data);
  std::vector<double>::const_iterator iter = squared.cbegin();
  std::advance(iter, 1);
  return std::reduce(iter, squared.cend(), 0.0);
};
const double
OrnsteinUhlenbeckLikelihoodComponentCalculator::calculateLagSumSquared(
    const std::vector<double>& data
) const {
  const std::vector<double> squared = valuesSquared(data);
  std::vector<double>::const_iterator iter = squared.cend();
  std::advance(iter, -1);
  return std::reduce(squared.cbegin(), iter, 0.0);
};
const double
OrnsteinUhlenbeckLikelihoodComponentCalculator::calculateLeadLagSumProduct(
    const std::vector<double>& data
) const {
  std::vector<double>::const_iterator iter_lead = data.cbegin();
  std::advance(iter_lead, 1);

  std::vector<double>::const_iterator iter_lag_start = data.cbegin();
  std::vector<double>::const_iterator iter_lag_end = data.cend();
  std::advance(iter_lag_end, -1);

  return std::inner_product(iter_lag_start, iter_lag_end, iter_lead, 0.0);
};
const double OrnsteinUhlenbeckLikelihoodComponentCalculator::updateLeadSum(
    const double& lead_sum, const double& new_observation
) const {
  return lead_sum + new_observation;
};
const double OrnsteinUhlenbeckLikelihoodComponentCalculator::updateLagSum(
    const double& lag_sum, const double& last_observation
) const {
  return lag_sum + last_observation;
};
const double
OrnsteinUhlenbeckLikelihoodComponentCalculator::updateLeadSumSquared(
    const double& lead_sum_squared, const double& new_observation
) const {
  return lead_sum_squared + std::pow(new_observation, 2);
};
const double
OrnsteinUhlenbeckLikelihoodComponentCalculator::updateLagSumSquared(
    const double& lag_sum_squared, const double& last_observation
) const {
  return lag_sum_squared + std::pow(last_observation, 2);
};
const double
OrnsteinUhlenbeckLikelihoodComponentCalculator::updateLeadLagSumProduct(
    const double& lead_lag_sum_product,
    const double& new_observation,
    const double& last_observation
) const {
  return lead_lag_sum_product + (last_observation * new_observation);
};
OrnsteinUhlenbeckLikelihoodComponents
OrnsteinUhlenbeckLikelihood::calculateComponents(
    const std::vector<double>& data
) const {
  return OrnsteinUhlenbeckLikelihoodComponents{
      component_calculator.calculateLeadSum(data),
      component_calculator.calculateLagSum(data),
      component_calculator.calculateLeadSumSquared(data),
      component_calculator.calculateLagSumSquared(data),
      component_calculator.calculateLeadLagSumProduct(data),
      static_cast<uint32_t>(data.size())
  };
};
const OrnsteinUhlenbeckLikelihoodComponents
OrnsteinUhlenbeckLikelihood::updateComponents(
    const OrnsteinUhlenbeckLikelihoodComponents& components,
    const double& new_observation,
    const double& last_observation
) const {
  return OrnsteinUhlenbeckLikelihoodComponents{
      component_calculator.updateLeadSum(components.lead_sum, new_observation),
      component_calculator.updateLagSum(components.lag_sum, last_observation),
      component_calculator.updateLeadSumSquared(
          components.lead_sum_squared, new_observation
      ),
      component_calculator.updateLagSumSquared(
          components.lag_sum_squared, last_observation
      ),
      component_calculator.updateLeadLagSumProduct(
          components.lead_lag_sum_product, new_observation, last_observation
      ),
      components.n_obs + 1
  };
}
const double OrnsteinUhlenbeckLikelihoodComponentCalculator::calculateMu(
    const OrnsteinUhlenbeckLikelihoodComponents& components
) const {
  return ((components.lead_sum * components.lag_sum_squared) -
          (components.lag_sum * components.lead_lag_sum_product)) /
         ((components.n_obs *
           (components.lag_sum_squared - components.lead_lag_sum_product)) -
          (std::pow(components.lag_sum, 2) -
           (components.lead_sum * components.lag_sum)));
};
const double OrnsteinUhlenbeckLikelihoodComponentCalculator::calculateAlpha(
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
const double OrnsteinUhlenbeckLikelihoodComponentCalculator::calculateSigma(
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
  sigma *= (1.0 / components.n_obs);
  sigma *= ((2 * exp_alpha) / (1 - std::pow(exp_alpha, 2)));
  return sigma;
};
const OrnsteinUhlenbeckParameters
OrnsteinUhlenbeckLikelihood::calculateParameters(
    const OrnsteinUhlenbeckLikelihoodComponents& components
) const {
  const double mu = component_calculator.calculateMu(components);
  const double alpha = component_calculator.calculateAlpha(components, mu);
  const double sigma =
      component_calculator.calculateSigma(components, mu, alpha);
  return OrnsteinUhlenbeckParameters{mu, alpha, sigma};
};
