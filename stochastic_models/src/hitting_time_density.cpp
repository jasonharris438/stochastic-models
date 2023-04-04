#include "numeric_utils/integration.h"
#include "sde/stochastic_model.h"
double integrate_hitting_time_density(double x, void* model) {
    StochasticModel* m = static_cast<StochasticModel*>(model);
    double value = m->HittingTimeDensityCore(x);
    return value;
}
const double hitting_time_density(double& x, model_func fn, void* model,
                                  double& first, double& second) {
    return integration(fn, model, second, x) /
           integration(fn, model, second, first);
}
