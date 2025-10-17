
#include "stochastic_models/likelihood/general_linear_likelihood.h"

#include "stochastic_models/numeric_utils/helpers.h"

#include <cmath>
#include <numeric>

GeneralLinearLikelihoodComponents::GeneralLinearLikelihoodComponents(
    const double lag_squared,
    const double lead_lag_inner_product,
    const double sigma_kernel,
    const uint32_t n_obs
)
    : lag_squared(lag_squared), lead_lag_inner_product(lead_lag_inner_product),
      sigma_kernel(sigma_kernel), n_obs(n_obs) {}

const double GeneralLinearLikelihood::calculateLeadLagInnerProduct(
    const std::vector<double>& data
) const {
  std::vector<double>::const_iterator iter_y = data.cbegin();
  std::advance(iter_y, 1);

  std::vector<double>::const_iterator iter_x_start = data.cbegin();
  std::vector<double>::const_iterator iter_x_end = data.cend();
  std::advance(iter_x_end, -1);

  return std::inner_product(iter_x_start, iter_x_end, iter_y, 0.0);
};
const double GeneralLinearLikelihood::calculateLagSquared(
    const std::vector<double>& data
) const {
  const std::vector<double> squared = valuesSquared(data);
  std::vector<double>::const_iterator iter = squared.cbegin();
  std::advance(iter, 1);
  return std::reduce(iter, squared.cend(), 0.0);
};
const double GeneralLinearLikelihood::calculateSigmaKernel(
    const u_int32_t& n_observations, const double& sigma
) const {
  return std::pow(sigma * n_observations, 2);
}
const double
GeneralLinearLikelihood::calculateMu(const std::vector<double>& data) const {
  return std::log(
      calculateLeadLagInnerProduct(data) / calculateLagSquared(data)
  );
};
const double GeneralLinearLikelihood::calculateSigma(
    const std::vector<double>& data, const double& mu
) const {
  const double exp_mu = std::exp(mu);
  double squared_diff_accumulated = 0.0;

  if (data.size() > 0) {
    for (size_t idx_lead = 1; idx_lead < data.size(); idx_lead++) {
      size_t idx_lag = idx_lead - 1;
      const double lag_factored = data.at(idx_lag) * exp_mu;

      squared_diff_accumulated += std::pow(data.at(idx_lead) - lag_factored, 2);
    }
    return std::sqrt(squared_diff_accumulated / data.size());
  } else {
    return squared_diff_accumulated;
  }
};
const double GeneralLinearLikelihood::calculateConditionalVariance(
    const double& sigma, const double& mu
) const {
  return (2 * sigma * mu) / (std::exp(2 * mu) - std::exp(mu));
};
const GeneralLinearParameters GeneralLinearLikelihood::calculateParameters(
    const std::vector<double>& data
) const {
  const double mu = calculateMu(data);
  const double sigma = calculateSigma(data, mu);
  return GeneralLinearParameters{mu, sigma};
}

const GeneralLinearLikelihoodComponents
GeneralLinearLikelihood::calculateComponents(
    const std::vector<double>& data, const GeneralLinearParameters& params
) const {
  const double lead_lag_inner_product = calculateLeadLagInnerProduct(data);
  const double lag_squared = calculateLagSquared(data);
  const uint32_t n_obs = static_cast<uint32_t>(data.size());
  const double sigma_kernel_value = calculateSigmaKernel(n_obs, params.sigma);
  return GeneralLinearLikelihoodComponents{
      lag_squared, lead_lag_inner_product, sigma_kernel_value, n_obs
  };
}
