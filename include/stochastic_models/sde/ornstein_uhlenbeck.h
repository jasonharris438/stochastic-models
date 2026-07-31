#ifndef STOCHASTIC_MODELS_SDE_ORNSTEIN_UHLENBECK_H
#define STOCHASTIC_MODELS_SDE_ORNSTEIN_UHLENBECK_H
#include "stochastic_models/sde/stochastic_model.h"

/**
 * @brief Handles fitting, evaluating, and simulating specifically the
 * Ornstein-Uhlenbeck model specification.
 */
class OrnsteinUhlenbeckModel : public StochasticModel::StochasticModel {
private:
  /**
   * @brief The mean of the Ornstein-Uhlenbeck model.
   *
   */
  const double mu;
  /**
   * @brief The mean-reverting velocity of the Ornstein-Uhlenbeck model.
   *
   */
  const double alpha;
  /**
   * @brief The standard deviation of the Ornstein-Uhlenbeck model.
   *
   */
  const double sigma;

public:
  OrnsteinUhlenbeckModel();
  OrnsteinUhlenbeckModel(
      const double mu, const double alpha, const double sigma
  );
  OrnsteinUhlenbeckModel(const OrnsteinUhlenbeckModel& other);
  /**
   * @brief Return a heap-allocated copy of the model (virtual constructor).
   * @return const OrnsteinUhlenbeckModel* Pointer to the new instance.
   */
  const OrnsteinUhlenbeckModel* clone() const override;
  /**
   * @brief Returns the model mean.
   *
   * @return const double The model mean.
   */
  const double getMean() const override;
  /**
   * @brief Returns the unconditional variance of the model.
   *
   * @return const double The model unconditional variance.
   */
  const double getUnconditionalVariance() const override;
  /**
   * @brief Conditional variance of the process over a step of length t:
   * sigma^2 * (1 - exp(-2 alpha t)) / (2 alpha), with alpha -> 0 limit
   * sigma^2 * t. See docs/derivations/sde-mle-derivations.md.
   *
   * @param t The step length.
   * @return const double The conditional variance over the step.
   */
  const double getConditionalVariance(const double t) const;
  /**
   * @brief Produces a simulation using the parameters mu, alpha, and sigma of
   * size provided in the method arguments. Uses coreEquation to produce the
   * model equation and the model distribution to produce .sample() values.
   *
   * @param start The value to start the simulation at.
   * @param size The number of values to simulate.
   * @return std::vector<double> A simulated model series.
   * @throws std::invalid_argument If `size < 1`.
   */
  std::vector<double> Simulate(
      const double start, const unsigned int& size, const unsigned int& t
  ) const override;
  /**
   * @brief Exact-scheme step of the Ornstein-Uhlenbeck process: conditional-
   * mean drift plus a Gaussian increment with the exact conditional standard
   * deviation. See docs/derivations/sde-mle-derivations.md.
   *
   * @param x The current value of the series.
   * @param noise The random Gaussian noise to add to the series.
   * @param t The time increment of a single step.
   * @return const double The next value in the series.
   */
  const double coreEquation(
      const double& x, const double& noise, const unsigned int& t
  ) const override;

  ~OrnsteinUhlenbeckModel() override;
};
#endif // STOCHASTIC_MODELS_SDE_ORNSTEIN_UHLENBECK_H
