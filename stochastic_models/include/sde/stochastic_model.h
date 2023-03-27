#ifndef _STOCHASTIC_MODEL_H_
#define _STOCHASTIC_MODEL_H_
#include "distributions/gaussian.h"
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
     * Simulates size many random draws from the CoreEquation implemented in the
     * child class. Noise is created from a standard normal distribution.
     *
     * @param start The starting point of the simulation.
     * @param size How many samples to draw. Defaults to 1.
     * @returns Random values drawn from CoreEquation.
     */
    virtual std::vector<double> Simulate(const double start,
                                         const unsigned int size) const = 0;

    /**
     * Implements the core model equation defined in the child class that
     * inherits from this base class.
     *
     * @param x The point at which to evaluate the core equation.
     * @param noise Gaussian noise to apply to the value produced.
     * @param t The size of the time steps (models are typically discretized).
     * @returns Core equation evaluated at x.
     */
    virtual const double CoreEquation(const double& x, const double noise,
                                      const unsigned int t = 1) const = 0;
    /**
     * @brief Implements the core components of the hitting time density
     * function for the model.
     *
     * @param x The point at which to evaluate the hitting time density function
     * core.
     * @return const double The hitting time density function core evaluated at
     * x.
     */
    virtual const double HittingTimeDensityCore(const double& x) const = 0;
    virtual ~StochasticModel() = 0;
};
#endif  // _STOCHASTIC_MODEL_H_
