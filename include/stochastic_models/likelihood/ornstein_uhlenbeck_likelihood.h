#ifndef STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_LIKELIHOOD_H
#define STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_LIKELIHOOD_H

/**
 * @file include/stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h
 * @brief Maximum likelihood helpers and calculators for the Ornstein-Uhlenbeck
 * model.
 * @ingroup likelihood
 * @details
 * This header exposes utilities to compute the summarising components used by
 * the Ornstein-Uhlenbeck maximum likelihood estimation routines and a small API
 * to compute the final parameters (mu, alpha, sigma) from those components. The
 * implementation is designed to allow efficient calculation of sums over the
 * lead/lag series and to support incremental updates.
 *
 * @note The functions operating on `OrnsteinUhlenbeckLikelihoodComponents`
 * expect `components.n_obs > 0` where indicated. Consumers should validate
 * input series length before calling the higher-level compute functions.
 */

#include <cstdint>
#include <vector>

struct OrnsteinUhlenbeckParameters {
  double mu;
  double alpha;
  double sigma;
};

struct OrnsteinUhlenbeckLikelihoodComponents {
  double lead_sum;
  double lag_sum;
  double lead_sum_squared;
  double lag_sum_squared;
  double lead_lag_sum_product;
  uint32_t n_obs;
};

class OrnsteinUhlenbeckLikelihoodComponentCalculator {
public:
  /**
   * @brief Calculates the sum of the lead series used in likelihood
   * calculations.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateLeadSum(const std::vector<double>& data) const;
  /**
   * @brief Calculates the sum of the lagged series used in likelihood
   * calculations.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateLagSum(const std::vector<double>& data) const;
  /**
   * @brief Calculates the sum of squares of the lead series used in likelihood
   * calculations.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateLeadSumSquared(const std::vector<double>& data) const;
  /**
   * @brief Calculates the sum of squares of the lagged series used in
   * likelihood calculations.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double calculateLagSumSquared(const std::vector<double>& data) const;
  /**
   * @brief Calculates the sum of the products of the lead and lagged series
   * used in likelihood calculations.
   *
   * @param data Data series used to calculate.
   * @return const double Likelihood component value.
   */
  const double
  calculateLeadLagSumProduct(const std::vector<double>& data) const;
  /**
   * @brief Incrementally update the sum of the lead series.
   *
   * @param lead_sum Current accumulated lead sum.
   * @param new_observation The new time t observation to add.
   * @return const double Updated lead sum.
   */
  const double
  updateLeadSum(const double& lead_sum, const double& new_observation) const;
  /**
   * @brief Incrementally update the sum of the lag series.
   *
   * @param lag_sum Current accumulated lag sum.
   * @param last_observation The new time t-1 observation to add.
   * @return const double Updated lag sum.
   */
  const double
  updateLagSum(const double& lag_sum, const double& last_observation) const;
  /**
   * @brief Incrementally update the sum of squares of the lead series.
   *
   * @param lead_sum_squared Current accumulated lead sum of squares.
   * @param new_observation The new time t observation to add.
   * @return const double Updated lead sum of squares.
   */
  const double updateLeadSumSquared(
      const double& lead_sum_squared, const double& new_observation
  ) const;
  /**
   * @brief Incrementally update the sum of squares of the lag series.
   *
   * @param lag_sum_squared Current accumulated lag sum of squares.
   * @param last_observation The new time t-1 observation to add.
   * @return const double Updated lag sum of squares.
   */
  const double updateLagSumSquared(
      const double& lag_sum_squared, const double& last_observation
  ) const;
  /**
   * @brief Incrementally update the sum of the products of the lead and lag
   * series.
   *
   * @param lead_lag_sum_product Current accumulated lead-lag sum product.
   * @param new_observation The new time t observation to add.
   * @param last_observation The new time t-1 observation to add.
   * @return const double Updated lead-lag sum product.
   */
  const double updateLeadLagSumProduct(
      const double& lead_lag_sum_product,
      const double& new_observation,
      const double& last_observation
  ) const;
  /**
   * @brief Calculates the maximum likelihood value for the model parameter
   * mu using the values in the components member.
   *
   * @param components Summary values of a data series used to calculate the
   * maximum likelihood parameters.
   * @return const double The value for mu parameter.
   */
  const double
  calculateMu(const OrnsteinUhlenbeckLikelihoodComponents& components) const;

  /**
   * @brief Calculates the maximum likelihood value for the model parameter
   * alpha using the values in the components member and the existing mu
   * value.
   *
   * @param components Summary values of a data series used to calculate the
   * maximum likelihood parameters.
   * @param mu The maximum likelihood parmater mu for the model.
   * @return const double The value for alpha parameter.
   */
  const double calculateAlpha(
      const OrnsteinUhlenbeckLikelihoodComponents& components, const double& mu
  ) const;
  /**
   * @brief Calculates the maximum likelihood value for the model parameter
   * sigma using the values in the components member and the existing mu and
   * alpha values.
   *
   * @param components Summary values of a data series used to calculate the
   * maximum likelihood parameters.
   * @param mu The maximum likelihood parmater mu for the model.
   * @param alpha The maximum likelihood parmater alpha for the model.
   * @return const double The value for sigma parameter.
   */
  const double calculateSigma(
      const OrnsteinUhlenbeckLikelihoodComponents& components,
      const double& mu,
      const double& alpha
  ) const;
};

/**
 * @brief Handles calculating maximum likelihood parameters specifically under
 * Ornstein-Uhlenbeck processes assumptions for a data series.
 *
 */
class OrnsteinUhlenbeckLikelihood {
private:
  /**
   * @brief Component calculator instance to use when calculating likelihood
   * components from a data series.
   */
  OrnsteinUhlenbeckLikelihoodComponentCalculator component_calculator;

public:
  /**
   * @brief Calculates the internal likelihood equation components from a data
   * series.
   *
   * @param data Data series used to calculate maximum likelihood parameters.
   * @return OrnsteinUhlenbeckLikelihoodComponents Struct with precomputed
   *         sums (lead/lag sums, squared sums, cross-products) and observation
   *         count `n_obs`.
   * @pre `data.size() >= 2` (requires at least one lead/lag pair).
   * @note This method is deterministic and performs no allocations beyond the
   *       returned struct.
   */
  OrnsteinUhlenbeckLikelihoodComponents
  calculateComponents(const std::vector<double>& data) const;

  /**
   * @brief Convert precomputed components into OU model parameters.
   *
   * Uses closed-form relationships derived from the log-likelihood to compute
   * the MLE for `mu`, `alpha`, and `sigma`.
   *
   * @param components Precomputed series summary values generated by
   *                   `calculateComponents`.
   * @return OrnsteinUhlenbeckParameters Structure containing {mu, alpha,
   * sigma}.
   * @pre `components.n_obs > 0`.
   * @post `result.sigma >= 0` unless inputs are degenerate.
   */
  const OrnsteinUhlenbeckParameters calculateParameters(
      const OrnsteinUhlenbeckLikelihoodComponents& components
  ) const;
  /**
   * @brief Incrementally update the components with a single new observation.
   *
   * This method allows online usage: given an existing components struct and a
   * (new_observation, last_observation) pair, it returns a new components
   * struct with updated sums and `n_obs`.
   *
   * @param components Current components state.
   * @param new_observation The new observed value X_t.
   * @param last_observation The previous observed value X_{t-1}.
   * @return OrnsteinUhlenbeckLikelihoodComponents Updated components with
   *         `n_obs` incremented by 1.
   * @pre components.n_obs >= 0
   */
  const OrnsteinUhlenbeckLikelihoodComponents updateComponents(
      const OrnsteinUhlenbeckLikelihoodComponents& components,
      const double& new_observation,
      const double& last_observation
  ) const;
};

#endif // STOCHASTIC_MODELS_LIKELIHOOD_ORNSTEIN_UHLENBECK_LIKELIHOOD_H
