#ifndef _TYPES_H_
#define _TYPES_H_
/**
 * @brief A function pointer to a function that takes a double and a pointer
 * to void. This type is used to ensure compatability with the GSL
 * integration API.
 *
 */
typedef double (*model_func)(double x, void* model);
#endif  // _TYPES_H_
