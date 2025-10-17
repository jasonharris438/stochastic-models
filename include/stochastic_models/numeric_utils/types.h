#ifndef STOCHASTIC_MODELS_NUMERIC_UTILS_TYPES_H
#define STOCHASTIC_MODELS_NUMERIC_UTILS_TYPES_H
/**
 * @brief A function pointer to a function that takes a double and a pointer
 * to void. This type is used to ensure compatability with the GSL
 * integration API.
 *
 */
typedef double (*ModelFunc)(double x, void* model);
#endif // STOCHASTIC_MODELS_NUMERIC_UTILS_TYPES_H
