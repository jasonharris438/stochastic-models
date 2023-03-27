#include <gsl/gsl_integration.h>
#include <math.h>
#include <stdio.h>

#include "numeric_utils/integration.h"
#include "sde/ornstein_uhlenbeck.h"
double f(double x, void* p) {
    struct params* parameters = (struct params*)p;
    double alpha = parameters->alpha;
    double mu = parameters->mu;
    double sigma = parameters->sigma;
    OrnsteinUhlenbeckModel model = OrnsteinUhlenbeckModel(mu, alpha, sigma);
    return model.HittingTimeDensityCore(x);
}
const double integration(double& alpha, double& mu, double& sigma,
                         double& lower, double& upper) {
    // https :  // www.gnu.org/software/gsl/doc/html/integration.html
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    double result, error;
    params p = {alpha, mu, sigma};

    gsl_function F;
    F.function = &f;
    F.params = &p;

    gsl_integration_qags(&F, lower, upper, 0, 1e-7, 1000, w, &result, &error);

    const double value = result;

    gsl_integration_workspace_free(w);

    return value;
}
const double hitting_time_density(double& x, double& alpha, double& mu,
                                  double& sigma, double& first,
                                  double& second) {
    return integration(alpha, mu, sigma, second, x) /
           integration(alpha, mu, sigma, second, first);
}