#ifndef _INTEGRATION_H_
#define _INTEGRATION_H_
#include "numeric_utils/types.h"
/**
 * @brief Integrates the function f over a given interval.
 *
 * @param fn A pointer to the function to integrate.
 * @param model A pointer to the model instance that is being used
 * to perform integration.
 * @param lower A reference to the lower bound of the integration interval.
 * @param upper A reference to the upper bound of the integration interval.
 * @return const double. The value of the integral of F.
 */
const double integration(model_func fn, void* model, double& lower,
                         double& upper);
#endif  // _INTEGRATION_H_
