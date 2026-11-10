
#include "stochastic_models/likelihood/general_linear_likelihood.h"

#include "stochastic_models/numeric_utils/helpers.h"

#include <cmath>
#include <iterator>
#include <numeric>

const double GeneralLinearLikelihoodComponentCalculator::calculateSeriesMean(
    const double& numerator, const double& denominator
) const {
  if (denominator == 0 || numerator == 0) {
    return 0.0;
  } else {
    return numerator / denominator;
  }
}

const double
GeneralLinearLikelihoodComponentCalculator::calculateLeadLagInnerProduct(
    const std::vector<double>& data
) const {
  check_minimum_observations(data, 2, "calculateLeadLagInnerProduct");
  std::vector<double>::const_iterator iter_y = data.cbegin();
  std::advance(iter_y, 1);

  std::vector<double>::const_iterator iter_x_start = data.cbegin();
  std::vector<double>::const_iterator iter_x_end = data.cend();
  std::advance(iter_x_end, -1);

  return std::inner_product(iter_x_start, iter_x_end, iter_y, 0.0);
};
const double GeneralLinearLikelihoodComponentCalculator::calculateLagSquared(
    const std::vector<double>& data
) const {
  check_minimum_observations(data, 2, "calculateLagSquared");
  const std::vector<double> squared = valuesSquared(data);
  return std::reduce(squared.cbegin(), std::prev(squared.cend()), 0.0);
};
const double GeneralLinearLikelihoodComponentCalculator::calculateSquaredError(
    const std::vector<double>& data, const double& mu
) const {
  check_minimum_observations(data, 2, "calculateSquaredError");
  const double exp_mean = std::exp(mu);
  double squared_diff_accumulated = 0.0;
  for (size_t idx_lead = 1; idx_lead < data.size(); idx_lead++) {
    size_t idx_lag = idx_lead - 1;
    const double lag_factored = data.at(idx_lag) * exp_mean;

    squared_diff_accumulated += std::pow(data.at(idx_lead) - lag_factored, 2);
  }
  return squared_diff_accumulated;
};
const double
GeneralLinearLikelihoodComponentCalculator::updateLeadLagInnerProduct(
    const double& lead_lag_inner_product,
    const double& new_observation,
    const double& last_observation
) const {
  return lead_lag_inner_product + (last_observation * new_observation);
}
const double GeneralLinearLikelihoodComponentCalculator::updateLagSquared(
    const double& lag_squared, const double& last_observation
) const {
  return lag_squared + std::pow(last_observation, 2);
}
const double GeneralLinearLikelihoodComponentCalculator::updateSquaredError(
    const double& squared_error,
    const double& new_observation,
    const double& last_observation,
    const double& prior_mean,
    const double& lag_squared_old,
    const double& lag_squared_new
) const {
  if (lag_squared_old == 0.0) {
    return squared_error;
  }
  const double residual = new_observation - prior_mean * last_observation;
  return squared_error +
         std::pow(residual, 2) * lag_squared_old / lag_squared_new;
}
const double GeneralLinearLikelihoodComponentCalculator::calculateMu(
    const GeneralLinearLikelihoodComponents& components
) const {
  return std::log(calculateSeriesMean(
      components.lead_lag_inner_product, components.lag_squared
  ));
};
const double GeneralLinearLikelihoodComponentCalculator::calculateSigma(
    const GeneralLinearLikelihoodComponents& components, const double& mu
) const {
  if (components.squared_error == 0.0 || components.n_obs < 2) {
    return 0.0;
  }
  const double n_pairs = components.n_obs - 1.0;
  // Transition variance to squared diffusion coefficient; the removable
  // singularity at mu == 0 has limit 1.
  const double conversion = (mu == 0.0) ? 1.0 : (2 * mu) / std::expm1(2 * mu);
  return std::sqrt(components.squared_error / n_pairs * conversion);
};
const GeneralLinearParameters GeneralLinearLikelihood::calculateParameters(
    const GeneralLinearLikelihoodComponents& components
) const {
  const double mu = component_calculator.calculateMu(components);
  const double sigma = component_calculator.calculateSigma(components, mu);
  return GeneralLinearParameters{mu, sigma};
}

const GeneralLinearLikelihoodComponents
GeneralLinearLikelihood::calculateComponents(
    const std::vector<double>& data
) const {
  check_minimum_observations(data, 2, "General linear maximum likelihood");
  const uint32_t n_obs = static_cast<uint32_t>(data.size());
  const double lead_lag_inner_product =
      component_calculator.calculateLeadLagInnerProduct(data);
  const double lag_squared = component_calculator.calculateLagSquared(data);
  return GeneralLinearLikelihoodComponents{
      lag_squared, lead_lag_inner_product,
      component_calculator.calculateSquaredError(
          data, std::log(component_calculator.calculateSeriesMean(
                    lead_lag_inner_product, lag_squared
                ))
      ),
      n_obs
  };
}
const GeneralLinearLikelihoodComponents
GeneralLinearLikelihood::updateComponents(
    const GeneralLinearLikelihoodComponents& components,
    const double& new_observation,
    const double& last_observation
) const {
  const double prior_mean = component_calculator.calculateSeriesMean(
      components.lead_lag_inner_product, components.lag_squared
  );
  const double lag_squared = component_calculator.updateLagSquared(
      components.lag_squared, last_observation
  );
  const double lead_lag_inner_product =
      component_calculator.updateLeadLagInnerProduct(
          components.lead_lag_inner_product, new_observation, last_observation
      );
  const double squared_error = component_calculator.updateSquaredError(
      components.squared_error, new_observation, last_observation, prior_mean,
      components.lag_squared, lag_squared
  );
  return GeneralLinearLikelihoodComponents{
      lag_squared, lead_lag_inner_product, squared_error, components.n_obs + 1
  };
}
