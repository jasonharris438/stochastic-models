#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/hitting_times/hitting_time_ornstein_uhlenbeck.h"
#include "stochastic_models/numeric_utils/helpers.h"
#include "stochastic_models/numeric_utils/integration.h"
double integrateHittingTimeDensity(double x, void* model) {
  HittingTimeOrnsteinUhlenbeck* m =
      static_cast<HittingTimeOrnsteinUhlenbeck*>(model);
  double value = m->hittingTimeDensityCore(x);
  return value;
}
const double hittingTimeDensity(
    double& x, ModelFunc fn, void* model, double& first, double& second
) {
  const double numerator = adaptiveIntegration(fn, model, second, x);
  const double denominator = adaptiveIntegration(fn, model, second, first);

  if (denominator == 0.0) {
    throw ZeroDivError("hittingTimeDensity denominator integral is zero.");
  }

  const double density = numerator / denominator;
  check_finite_result(density, "hittingTimeDensity");

  return density;
}
