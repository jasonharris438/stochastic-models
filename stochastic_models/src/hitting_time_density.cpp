#include <gsl/gsl_integration.h>
#include <math.h>
#include <stdio.h>

#include "numeric_utils/integration.h"
#include "sde/ornstein_uhlenbeck.h"
double f(double x, void *p) {
    struct params *parameters = (struct params *)p;
    double alpha = parameters->alpha;
    double mu = parameters->mu;
    double sigma = parameters->sigma;
    OrnsteinUhlenbeckModel model = OrnsteinUhlenbeckModel(mu, alpha, sigma);
    return model.HittingTimeDensityCore(x);
}
double integration() {
    // https :  // www.gnu.org/software/gsl/doc/html/integration.html
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    double result, error;
    double expected = -4.0;
    double alpha = 5.1;
    double mu = 0.996;
    double sigma = 1.1;
    double upper = 1.05;
    double lower = 0.8;
    params p = {alpha, mu, sigma};

    gsl_function F;
    F.function = &f;
    F.params = &p;

    gsl_integration_qags(&F, lower, upper, 0, 1e-7, 1000, w, &result, &error);

    printf("result          = % .18f\n", result);
    printf("exact result    = % .18f\n", expected);
    printf("estimated error = % .18f\n", error);
    printf("actual error    = % .18f\n", result - expected);
    printf("intervals       = %zu\n", w->size);

    double value = result;

    gsl_integration_workspace_free(w);

    return value;
}
