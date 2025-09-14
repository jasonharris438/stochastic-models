#ifndef ORNSTEIN_UHLENBECK_H
#define ORNSTEIN_UHLENBECK_H
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
    OrnsteinUhlenbeckModel(const double mu, const double alpha,
                           const double sigma);
    OrnsteinUhlenbeckModel(const OrnsteinUhlenbeckModel& other);
    /**
     * @brief Construct a new OrnsteinUhlenbeckModel object and return on heap
     * memory using the class' copy constructor in the caller instance
     *
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
     * @brief Produces a simulation using the parameters mu, alpha, and sigma of
     * size provided in the method arguments. Uses coreEquation to produce the
     * model equation and the model distribution to produce .sample() values.
     *
     * @param start The value to start the simulation at.
     * @param size The number of values to simulate.
     * @return std::vector<double> A simulated model series.
     */
    std::vector<double> Simulate(const double start, const unsigned int& size,
                                 const unsigned int& t) const override;
    /**
     * @brief Uses the Euler–Maruyama method for the approximate numerical
     * solution of the Ornstein-Uhlenbeck process.
     *
     * @param x The current value of the series.
     * @param noise The random Gaussian noise to add to the series.
     * @param t The time increment of a single step.
     * @return const double The next value in the series.
     */
    const double coreEquation(const double& x, const double& noise,
                              const unsigned int& t) const override;

    ~OrnsteinUhlenbeckModel() override;
};
#endif  // ORNSTEIN_UHLENBECK_H
