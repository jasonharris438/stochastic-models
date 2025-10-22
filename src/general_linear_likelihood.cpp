
#include "stochastic_models/likelihood/general_linear_likelihood.h"

#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/numeric_utils/helpers.h"

#include <cmath>
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
  const std::vector<double> squared = valuesSquared(data);
  std::vector<double>::const_iterator iter = squared.cbegin();
  std::advance(iter, 1);
  return std::reduce(iter, squared.cend(), 0.0);
};
const double GeneralLinearLikelihoodComponentCalculator::calculateSquaredError(
    const std::vector<double>& data, const double& mu
) const {
  const double exp_mean = std::exp(mu);
  double squared_diff_accumulated = 0.0;
  if (data.size() > 1) {
    for (size_t idx_lead = 1; idx_lead < data.size(); idx_lead++) {
      size_t idx_lag = idx_lead - 1;
      const double lag_factored = data.at(idx_lag) * exp_mean;

      squared_diff_accumulated += std::pow(data.at(idx_lead) - lag_factored, 2);
    }
  } else {
    throw InvalidNumberObservationsError(
        "Number of observations must be greater than 1."
    );
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
    const double& mean,
    const uint32_t& n_observations
) const {
  const double n_observations_new = n_observations + 1.0;
  const double n_observations_ratio = n_observations / n_observations_new;
  return squared_error +
         (n_observations_ratio *
          std::pow(new_observation - mean * last_observation, 2));
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
  const double sigma =
      (components.squared_error != 0.0 && components.n_obs > 0)
          ? std::sqrt(components.squared_error / components.n_obs)
          : 0.0;
  return sigma;
};
const double
GeneralLinearLikelihoodComponentCalculator::calculateConditionalVariance(
    const GeneralLinearParameters& parameters
) const {
  return (2 * parameters.sigma * parameters.mu) /
         (std::exp(2 * parameters.mu) - std::exp(parameters.mu));
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
const double GeneralLinearLikelihood::calculateConditionalVariance(
    const GeneralLinearParameters& parameters
) const {
  return component_calculator.calculateConditionalVariance(parameters);
}
const GeneralLinearLikelihoodComponents
GeneralLinearLikelihood::updateComponents(
    const GeneralLinearLikelihoodComponents& components,
    const double& new_observation,
    const double& last_observation
) const {
  const double lag_squared = component_calculator.updateLagSquared(
      components.lag_squared, last_observation
  );
  const double lead_lag_inner_product =
      component_calculator.updateLeadLagInnerProduct(
          components.lead_lag_inner_product, new_observation, last_observation
      );
  const double mean = component_calculator.calculateSeriesMean(
      lead_lag_inner_product, lag_squared
  );
  const double squared_error = component_calculator.updateSquaredError(
      components.squared_error, new_observation, last_observation, mean,
      components.n_obs
  );
  return GeneralLinearLikelihoodComponents{
      lag_squared, lead_lag_inner_product, squared_error, components.n_obs + 1
  };
}
