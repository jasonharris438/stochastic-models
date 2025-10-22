#ifndef STOCHASTIC_MODELS_ENTRYPOINTS_OU_MODEL_H
#define STOCHASTIC_MODELS_ENTRYPOINTS_OU_MODEL_H
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

/**
 * @brief Simulates an Ornstein-Uhlenbeck process of size provided in the
 * method arguments using the parameters mu, alpha, and sigma of the class
 * instance.
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param start The value to start the simulation at.
 * @param size The number of values to simulate.
 * @param t The time period to simulate over.
 * @return std::vector<double> A simulated Ornstein-Uhlenbeck process series.
 */
const std::vector<double> simulateOrnsteinUhlenbeck(
    const double& mu,
    const double& alpha,
    const double& sigma,
    const double start,
    const unsigned int& size,
    const unsigned int& t
);

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
const double hittingTimeDensityOrnsteinUhlenbeck(
    double x,
    const double mu,
    const double alpha,
    const double sigma,
    double first,
    double second
);
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

const std::vector<double> updateOuModel(
    const double mu,
    const double alpha,
    const double sigma,
    const double lead_sum,
    const double lag_sum,
    const double lead_sum_squared,
    const double lag_sum_squared,
    const double lead_lag_sum_product,
    const uint32_t n_obs,
    const double new_observation,
    const double last_observation

);

#endif // STOCHASTIC_MODELS_ENTRYPOINTS_OU_MODEL_H
