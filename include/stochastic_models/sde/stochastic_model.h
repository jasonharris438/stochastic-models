#ifndef STOCHASTIC_MODELS_SDE_STOCHASTIC_MODEL_H
#define STOCHASTIC_MODELS_SDE_STOCHASTIC_MODEL_H
#include "stochastic_models/distributions/gaussian.h"
/**
 * Stochastic Model base class that handles functionality for fitting,
 * analysing, and simulating statistical models. Should be treated as an
 * abstract class as most important functionality specific to a given model is
 * implemented in child classes.
 */
class StochasticModel {
protected:
  const GaussianDistribution* dist;

public:
  /**
   * @brief Construct a new StochasticModel object and return on heap memory
   * using the class' copy constructor in the caller instance
   *
   * @return const StochasticModel* Pointer to the new instance.
   */
  virtual const StochasticModel* clone() const = 0;
  /**
   * @brief Abstract method template for returning the model mean.
   *
   * @return const double The model mean.
   */
  virtual const double getMean() const = 0;
  /**
   * @brief Abstract method template for returning the unconditional variance
   * of the model.
   *
   * @return const double The model unconditional variance.
   */
  virtual const double getUnconditionalVariance() const = 0;
  /**
   * Simulates size many random draws from the coreEquation implemented in the
   * child class. Noise is created from a standard normal distribution.
   *
   * @param start The starting point of the simulation.
   * @param size How many samples to draw.
   * @param t The size of the time steps (models are typically discretized).
   * @returns Random values drawn from coreEquation.
   */
  virtual std::vector<double> Simulate(
      const double start, const unsigned int& size, const unsigned int& t
  ) const = 0;

  /**
   * Implements the core model equation defined in the child class that
   * inherits from this base class.
   *
   * @param x The point at which to evaluate the core equation.
   * @param noise Gaussian noise to apply to the value produced.
   * @param t The size of the time steps (models are typically discretized).
   * @returns Core equation evaluated at x.
   */
  virtual const double coreEquation(
      const double& x, const double& noise, const unsigned int& t
  ) const = 0;

  virtual ~StochasticModel() = 0;
};
#endif // STOCHASTIC_MODELS_SDE_STOCHASTIC_MODEL_H
