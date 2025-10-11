#include "stochastic_models/hitting_times/hitting_time_ornstein_uhlenbeck.h"
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
  return adaptiveIntegration(fn, model, second, x) /
         adaptiveIntegration(fn, model, second, first);
}
