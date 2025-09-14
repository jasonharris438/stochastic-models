#ifndef STOCHASTIC_models_H
#define STOCHASTIC_models_H

#include <string>
#include <unordered_map>
#include <vector>

#include "stochastic_models/kalman_filter/filter.h"
/**
 * @brief Calculate the optimal trading level exit value when a stop loss is
 * provided.
 *
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param stop_loss The stop loss value.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level exit value.
 */
const double optimalExitLevel(const double mu, const double alpha,
                              const double sigma, const double stop_loss,
                              const double r, const double c);
/**
 * @brief Calculate the optimal trading level exit value.
 *
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level exit value.
 */
const double optimalExitLevel(const double mu, const double alpha,
                              const double sigma, const double r,
                              const double c);
/**
 * @brief Calculate the optimal trading level exit value for an exponential
 * stochatic model.
 *
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level exit value.
 */
const double optimalExitLevelExponential(const double mu, const double alpha,
                                         const double sigma, const double r,
                                         const double c);
/**
 * @brief Calculate the lower bound of the optimal trading entry interval for an
 * exponential stochastic model.
 *
 * @param d_star The optimal trading level upper entry value.
 * @param b_star The optimal trading level exit value.
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level entry value.
 */
const double optimalEntryLevelLowerExponential(
    const double d_star, const double b_star, const double mu,
    const double alpha, const double sigma, const double r, const double c);
/**
 * @brief Calculate the lower bound of the optimal trading entry interval when a
 * stop loss is provided.
 *
 * @param d_star The optimal trading level upper entry value.
 * @param b_star The optimal trading level exit value.
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param stop_loss The stop loss value.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level entry value.
 */
const double optimalEntryLevelLower(const double d_star, const double b_star,
                                    const double mu, const double alpha,
                                    const double sigma, const double stop_loss,
                                    const double r, const double c);
/**
 * @brief Calculate the optimal trading level entry value for an exponential
 * stochatsic model.
 *
 * @param b_star The optimal trading level exit value.
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level entry value.
 */
const double optimalEntryLevelExponential(const double b_star, const double mu,
                                          const double alpha,
                                          const double sigma, const double r,
                                          const double c);
/**
 * @brief Calculate the optimal trading level entry value when a stop loss is
 * provided.
 *
 * @param b_star The optimal trading level exit value.
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param stop_loss The stop loss value.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level entry value.
 */
const double optimalEntryLevel(const double b_star, const double mu,
                               const double alpha, const double sigma,
                               const double stop_loss, const double r,
                               const double c);
/**
 * @brief Calculate the optimal trading level entry value.
 *
 * @param b_star The optimal trading level exit value.
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level entry value.
 */
const double optimalEntryLevel(const double b_star, const double mu,
                               const double alpha, const double sigma,
                               const double r, const double c);
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

/**
 * @brief Returns a JSON string containing an intialized kinetic components
 * state given the parameters provided.
 *
 * @param data_series The data series to initialize the filter with.
 * @param h A value determining the first and second derivative values of
 * the KCA system.
 * @param q A value determining the transition covariance of the KCA system.
 * @param system_dimensions A type containing the dimensions of the system
 * components.
 * @return const std::string The JSON string containing the initialized KCA
 * state.
 */
const std::string getInitializedKcaState(const std::vector<double> data_series,
                                         const double h, const double q,
                                         const std::string system_dimensions);
/**
 * @brief Takes a current kinetic components state and returns a JSON string
 * containing a single step update of the state and the value of the observation
 * provided.
 *
 * @param state A JSON string containing current state of the KCA system.
 * @param system_dimensions A type containing the dimensions of the system
 * components.
 * @param observation The next observation to update the system with.
 * @param innovation_sigma The sigma value of the innovation of the observed
 * data.
 * @return const std::string The JSON string containing the updated KCA state.
 */
const std::string getUpdatedKcaState(const std::string state,
                                     const std::string system_dimensions,
                                     const double observation,
                                     const double innovation_sigma);

#endif  // STOCHASTIC_models_H
