#include "numeric_utils/integration.h"

#include <gsl/gsl_integration.h>
const double integration(model_func fn, void* model, double& lower,
                         double& upper) {
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);

    double result, error;

    gsl_function F;
    F.function = *fn;
    F.params = model;

    gsl_integration_qags(&F, lower, upper, 0, 1e-7, 1000, w, &result, &error);

    const double value = result;

    gsl_integration_workspace_free(w);

    return value;
}
