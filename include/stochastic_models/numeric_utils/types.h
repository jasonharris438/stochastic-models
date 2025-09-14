#ifndef TYPES_H
#define TYPES_H
/**
 * @brief A function pointer to a function that takes a double and a pointer
 * to void. This type is used to ensure compatability with the GSL
 * integration API.
 *
 */
typedef double (*ModelFunc)(double x, void* model);
#endif  // TYPES_H
