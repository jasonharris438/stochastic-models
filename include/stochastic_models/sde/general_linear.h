#ifndef GENERAL_LINEAR_H
#define GENERAL_LINEAR_H
#include "stochastic_models/sde/stochastic_model.h"
/**
 * @brief General Linear Model class that handles fitting, evaluating, and
 * simulating specifically the General Linear SDE model specification.
 */
class GeneralLinearModel : public StochasticModel::StochasticModel {
private:
  /**
   * @brief The mean (mu parameter) of the General Linear model.
   *
   */
  const double mu;
  /**
   * @brief The standard deviation (sigma parameter) of the General Linear
   * model.
   *
   */
  const double sigma;

public:
  GeneralLinearModel();
  GeneralLinearModel(const double mu, const double sigma);
  GeneralLinearModel(const GeneralLinearModel& other);

  /**
   * @brief Construct a new GeneralLinearModel object and return on heap
   * memory using the class' copy constructor in the caller instance
   *
   * @return const GeneralLinearModel* Pointer to the new instance.
   */
  const GeneralLinearModel* clone() const override;

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
   * @brief Returns the conditional variance of the model.
   *
   * @return const double The model conditional variance.
   */
  const double getConditionalVariance() const;

  /**
   * @brief Produces a simulation using the parameters mu and sigma of
   * size provided in the method arguments. Uses coreEquation to produce the
   * model equation and the model distribution to produce .sample() values.
   *
   * @param start The value to start the simulation at.
   * @param size The number of values to simulate.
   * @return std::vector<double> A simulated model series.
   */
  std::vector<double> Simulate(
      const double start, const unsigned int& size, const unsigned int& t
  ) const override;
  /**
   * @brief Uses the Euler–Maruyama method for the approximate numerical
   * solution of the general linear SDE process.
   *
   * @param x The current value of the series.
   * @param noise The random Gaussian noise to add to the series.
   * @param t The time increment of a single step.
   * @return const double The next value in the series.
   */
  const double coreEquation(
      const double& x, const double& noise, const unsigned int& t
  ) const override;
  ~GeneralLinearModel() override;
};
#endif // GENERAL_LINEAR_H
