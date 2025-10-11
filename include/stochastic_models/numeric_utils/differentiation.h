#ifndef DIFFERENTIATION_H
#define DIFFERENTIATION_H
#include "stochastic_models/numeric_utils/types.h"
/**
 * @brief Computes the numerical derivative of the function @f[ F(x;r) @f] at
 * the point @f[ x @f] using an adaptive central difference algorithm with a
 * step-size of @f[ h @f]
 *
 * @param fn A pointer to the function to centrally difference.
 * @param model A pointer to the model instance that is being used
 * to perform central differencing.
 * @param x The point at which to perform differentiation.
 * @return const double. The value of the integral of @f[ F(x;r) @f].
 */
const double
adaptiveCentralDifferentiation(ModelFunc fn, void* model, double& x);
#endif // DIFFERENTIATION_H
