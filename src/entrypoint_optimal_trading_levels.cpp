/**
 * @file entrypoint_optimal_trading_levels.cpp
 * @brief Entry point module for optimal entry and exit trading levels
 * calculations.
 *
 */
#include <iostream>

#include "stochastic_models/entrypoints/optimal_trading_levels.h"
#include "stochastic_models/trading/exponential_mean_reversion.h"
#include "stochastic_models/trading/optimal_mean_reversion.h"
#include "stochastic_models/trading/trading_levels.h"
#include "stochastic_models/trading/trading_levels_exponential.h"

const double optimalExitLevel(const double mu, const double alpha,
                              const double sigma, const double stop_loss,
                              const double r, const double c) {
    // Object with outer/controller layer logic for trading levels.
    OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

    double value{0.0};
    try {
        // Calculate b*.
        value = tradingLevels.optimalExit(stop_loss, r, c);
    } catch (const std::runtime_error &e) {
        std::cout << "Exception " << e.what()
                  << " caught in optimalExitLevel with stop loss." << std::endl;
        throw;
    }
    return value;
}
const double optimalExitLevel(const double mu, const double alpha,
                              const double sigma, const double r,
                              const double c) {
    // Object with outer/controller layer logic for trading levels.
    OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

    double value{0.0};
    try {
        // Calculate b*.
        value = tradingLevels.optimalExit(r, c);
    } catch (const std::runtime_error &e) {
        std::cout << "Exception " << e.what()
                  << " caught in optimalExitLevel without stop loss."
                  << std::endl;
        throw;
    }
    return value;
}
const double optimalExitLevelExponential(const double mu, const double alpha,
                                         const double sigma,
                                         const double stop_loss, const double r,
                                         const double c) {
    // Object with outer/controller layer logic for trading levels.
    OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

    double value{0.0};
    try {
        // Calculate b*.
        value = tradingLevels.optimalExit(stop_loss, r, c);
    } catch (const std::runtime_error &e) {
        std::cout << "Exception " << e.what()
                  << " caught in optimalExitLevelExponential with stop loss."
                  << std::endl;
        throw;
    }
    return value;
}
const double optimalExitLevelExponential(const double mu, const double alpha,
                                         const double sigma, const double r,
                                         const double c) {
    // Object with outer/controller layer logic for trading levels.
    OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

    double value{0.0};
    try {
        // Calculate b*.
        value = tradingLevels.optimalExit(r, c);
    } catch (const std::runtime_error &e) {
        std::cout << "Exception " << e.what()
                  << " caught in optimalExitLevelExponential without stop loss."
                  << std::endl;
        throw;
    }
    return value;
}
const double optimalEntryLevel(const double b_star, const double mu,
                               const double alpha, const double sigma,
                               const double r, const double c) {
    // Object with outer/controller layer logic for trading levels.
    OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

    double value{0.0};
    try {
        // Calculate d*.
        value = tradingLevels.optimalEntry(b_star, r, c);
    } catch (const std::runtime_error &e) {
        std::cout << "Exception " << e.what()
                  << " caught in optimalEntryLevel without stop loss."
                  << std::endl;
        throw;
    }
    return value;
}
const double optimalEntryLevelExponential(const double b_star, const double mu,
                                          const double alpha,
                                          const double sigma, const double r,
                                          const double c) {
    // Object with outer/controller layer logic for trading levels.
    OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

    double value{0.0};
    try {
        // Calculate d*.
        value = tradingLevels.optimalEntry(b_star, r, c);
    } catch (const std::runtime_error &e) {
        std::cout
            << "Exception " << e.what()
            << " caught in optimalEntryLevelExponential without stop loss."
            << std::endl;
        throw;
    }
    return value;
}
const double optimalEntryLevelLowerExponential(
    const double d_star, const double b_star, const double mu,
    const double alpha, const double sigma, const double r, const double c) {
    // Object with outer/controller layer logic for trading levels.
    OrnsteinUhlenbeckTradingLevelsExponential tradingLevels(mu, alpha, sigma);

    double value{0.0};
    try {
        // Calculate d*.
        value = tradingLevels.optimalEntryLower(d_star, b_star, r, c);
    } catch (const std::runtime_error &e) {
        std::cout << "Exception " << e.what()
                  << " caught in optimalEntryLevelLowerExponential."
                  << std::endl;
        throw;
    }
    return value;
}
const double optimalEntryLevelLower(const double d_star, const double b_star,
                                    const double mu, const double alpha,
                                    const double sigma, const double stop_loss,
                                    const double r, const double c) {
    // Object with outer/controller layer logic for trading levels.
    OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

    double value{0.0};
    try {
        // Calculate d*.
        value =
            tradingLevels.optimalEntryLower(d_star, b_star, stop_loss, r, c);

    } catch (const std::runtime_error &e) {
        std::cout << "Exception " << e.what()
                  << " caught in optimalEntryLevelLower with stop loss."
                  << std::endl;
        throw;
    }
    return value;
}
const double optimalEntryLevel(const double b_star, const double mu,
                               const double alpha, const double sigma,
                               const double stop_loss, const double r,
                               const double c) {
    // Object with outer/controller layer logic for trading levels.
    OrnsteinUhlenbeckTradingLevels tradingLevels(mu, alpha, sigma);

    double value{0.0};
    try {
        // Calculate d*.
        value = tradingLevels.optimalEntry(b_star, stop_loss, r, c);
    } catch (const std::runtime_error &e) {
        std::cout << "Exception " << e.what()
                  << " caught in optimalEntryLevel with stop loss."
                  << std::endl;
        throw;
    }
    return value;
}
