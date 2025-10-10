#ifndef OU_MODEL_H
#define OU_MODEL_H
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @brief Calculates the probability of hitting the level first before second
 * given the parameters provided.
 *
 * @param x The starting point of the hitting time density calculation.
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param first The first level that is hit.
 * @param second The level first is assumed to be hit before.
 * @return const double The probability of hitting the level first before
 * second.
 */
const double hittingTimeDensityOrnsteinUhlenbeck(double x, const double mu,
                                                 const double alpha,
                                                 const double sigma,
                                                 double first, double second);
/**
 * @brief Calculate the maximum likelihood estimates of the
 * Ornstein-Uhlenbeck model parameters using the series in vec.
 *
 * @param vec The series to calculate the maximum likelihood estimates of
 * the Ornstein-Uhlenbeck model with.
 * @return const std::unordered_map<std::string, const double> The maximum
 * likelihood estimates of the Ornstein-Uhlenbeck model parameters.
 */
const std::unordered_map<std::string, const double>
ornsteinUhlenbeckMaximumLikelihood(const std::vector<double> vec);

#endif  // OU_MODEL_H
