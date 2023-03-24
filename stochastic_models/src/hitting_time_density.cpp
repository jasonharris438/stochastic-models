#include <gsl/gsl_integration.h>
#include <math.h>
#include <stdio.h>
/// @brief Parameters to be passed to the integrand function.
struct params {
    double alpha;
    double mu;
    double sigma;
};
/// The solution to the backward equation, the scale function @f[ S(x) @f]
/// @param x Point at which to evaluate the scale function.
/// @param p The parameters used in the scale function.
/// @return The value of the scale function at @f[ x @f].
double f(double x, void *p) {
    struct params *parameters = (struct params *)p;
    double alpha = parameters->alpha;
    double mu = parameters->mu;
    double sigma = parameters->sigma;
    return exp(x * alpha * (x - 2 * mu) / (pow(sigma, 2)));
}
// https :  // www.gnu.org/software/gsl/doc/html/integration.html
void integration() {
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

    gsl_integration_workspace_free(w);
}
