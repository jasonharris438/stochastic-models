#ifndef _ORNSTEIN_UHLENBECK_H_
#define _ORNSTEIN_UHLENBECK_H_
#include "stochastic_model.h"
/**
 * @brief Handles fitting, analysing, and simulating specifically the
 * Ornstein-Uhlenbeck model specification.
 */
class OrnsteinUhlenbeckModel : public StochasticModel {
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
    OrnsteinUhlenbeckModel(const double mu, const double alpha,
                           const double sigma);
    /**
     * @brief Produces a simulation using the parameters mu, alpha, and sigma of
     * size provided in the method arguments. Uses CoreEquation to produce the
     * model equation and the model distribution to produce .Sample() values.
     *
     * @param start The value to start the simulation at.
     * @param size The number of values to simulate.
     * @return std::vector<double> A simulated model series.
     */
    std::vector<double> Simulate(const double start,
                                 const unsigned int size) const override;
    /**
     * @brief Uses the Euler–Maruyama method for the approximate numerical
     * solution of the Ornstein-Uhlenbeck process.
     *
     * @param x The current value of the series.
     * @param noise The random Gaussian noise to add to the series.
     * @param t The time increment of a single step.
     * @return const double The next value in the series.
     */
    const double CoreEquation(const double x, const double noise,
                              const unsigned int t = 1) const override;
    /**
     * @brief Implements the core components of the hitting time density
     * function for the Ornstein-Uhlenbeck model.
     *
     * @param x The point at which to evaluate the hitting time density function
     * core.
     * @return const double The hitting time density function core evaluated at
     * x.
     */
    const double HittingTimeDensityCore(const double x) const override;
    ~OrnsteinUhlenbeckModel();
};
#endif  // _ORNSTEIN_UHLENBECK_H_
