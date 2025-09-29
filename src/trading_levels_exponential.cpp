#include "stochastic_models/trading/trading_levels_exponential.h"

#include <algorithm>
#include <iostream>

#include "stochastic_models/numeric_utils/helpers.h"
#include "stochastic_models/numeric_utils/solvers.h"
#include "stochastic_models/trading/trading_levels_params.h"

OrnsteinUhlenbeckTradingLevelsExponential::
    OrnsteinUhlenbeckTradingLevelsExponential(const double mu,
                                              const double alpha,
                                              const double sigma)
    : optimizer(std::make_unique<ExponentialMeanReversion>()),
      model(std::make_unique<OrnsteinUhlenbeckModel>(mu, alpha, sigma)),
      hitting_time_kernel(
          std::make_unique<HittingTimeOrnsteinUhlenbeck>(mu, alpha, sigma)) {}
const ExponentialMeanReversion*
OrnsteinUhlenbeckTradingLevelsExponential::getOptimizer() const {
    return optimizer.get();
};
const StochasticModel* OrnsteinUhlenbeckTradingLevelsExponential::getModel()
    const {
    return model.get();
};
const HittingTimeOrnsteinUhlenbeck*
OrnsteinUhlenbeckTradingLevelsExponential::getHittingTimeKernel() const {
    return hitting_time_kernel.get();
}
const ExponentialMeanReversion*
OrnsteinUhlenbeckTradingLevelsExponential::newOptimizer() const {
    return optimizer.get()->clone();
};
const StochasticModel* OrnsteinUhlenbeckTradingLevelsExponential::newModel()
    const {
    return model.get()->clone();
};
const HittingTimeOrnsteinUhlenbeck*
OrnsteinUhlenbeckTradingLevelsExponential::newHittingTimeKernel() const {
    return hitting_time_kernel.get()->clone();
}
const double OrnsteinUhlenbeckTradingLevelsExponential::optimalExitLowerBound(
    const double& r, const double& c) const {
    return std::max(getOptimizer()->L_star(getHittingTimeKernel(), r, c), c);
}
const double OrnsteinUhlenbeckTradingLevelsExponential::optimalExitUpperBound()
    const {
    return upperSolverBound(getModel());
}
const double OrnsteinUhlenbeckTradingLevelsExponential::optimalEntryLowerBound()
    const {
    return lowerSolverBound(getModel());
}
const double OrnsteinUhlenbeckTradingLevelsExponential::optimalExit(
    const double& stop_loss, const double& r, const double& c) const {
    // We need deep copies of the model and optimizer pointers to initialise the
    // params instance. This is because GSL requires a pointer to void and we
    // cannot use smart pointers with much benefit here. So we create deep
    // copies and then free that memory in the destructor of the
    // ExitLevelStopLossParams struct.
    void* params = new ExitLevelStopLossParams{getOptimizer()->clone(),
                                               getHittingTimeKernel()->clone(),
                                               stop_loss, r, c};
    ModelFunc fn = funcOptimalMeanReversionB;
    double value{0.0};
    try {
        double upper = optimalExitUpperBound();
        double lower = optimalExitLowerBound(r, c);
        value = brentSolver(fn, params, lower, upper);
    } catch (const std::exception& e) {
        std::cout << "Exception " << e.what()
                  << " caught in "
                     "OrnsteinUhlenbeckTradingLevelsExponential::optimalExit "
                     "with stop loss."
                  << std::endl;
        delete static_cast<ExitLevelStopLossParams*>(params);
        params = nullptr;
        throw;
    }
    // Then cast the void pointer back to the original type and free the memory.
    delete static_cast<ExitLevelStopLossParams*>(params);
    params = nullptr;
    return value;
}
const double OrnsteinUhlenbeckTradingLevelsExponential::optimalExit(
    const double& r, const double& c) const {
    // We need deep copies of the model and optimizer pointers to initialise the
    // params instance. This is because GSL requires a pointer to void and we
    // cannot use smart pointers with much benefit here. So we create deep
    // copies and then free that memory in the destructor of the ExitLevelParams
    // struct.
    void* params = new ExitLevelParams{optimizer->clone(),
                                       hitting_time_kernel->clone(), r, c};
    ModelFunc fn = funcOptimalMeanReversionB;

    double value{0.0};
    try {
        double upper = optimalExitUpperBound();
        double lower = optimalExitLowerBound(r, c);
        value = brentSolver(fn, params, lower, upper);
    } catch (const std::exception& e) {
        std::cout << "Exception " << e.what()
                  << " caught in "
                     "OrnsteinUhlenbeckTradingLevelsExponential::optimalExit."
                  << std::endl;
        delete static_cast<ExitLevelParams*>(params);
        params = nullptr;
        throw;
    }
    // Then cast the void pointer back to the original type and free the memory.
    delete static_cast<ExitLevelParams*>(params);
    params = nullptr;
    return value;
}
const double OrnsteinUhlenbeckTradingLevelsExponential::optimalEntryLower(
    const double& d_star, const double& b_star, const double& r,
    const double& c) const {
    // We need deep copies of the model and optimizer pointers to initialise the
    // params instance. This is because GSL requires a pointer to void and we
    // cannot use smart pointers with much benefit here. So we create deep
    // copies and then free that memory in the destructor of the
    // EntryLevelParams struct.
    void* params = new EntryLevelParams{newOptimizer(), newHittingTimeKernel(),
                                        b_star, r, c};
    ModelFunc fn = funcOptimalMeanReversionA;

    double value{0.0};
    try {
        double lower = optimalEntryLowerBound();
        double upper = d_star;
        value = brentSolver(fn, params, lower, upper);
    } catch (const std::exception& e) {
        std::cout
            << "Exception " << e.what()
            << " caught in "
               "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntryLower "
               "without stop loss."
            << std::endl;
        delete static_cast<ExitLevelParams*>(params);
        params = nullptr;
        throw;
    }

    // Then cast the void pointer back to the original type and free the memory.
    delete static_cast<EntryLevelParams*>(params);
    params = nullptr;
    return value;
}
const double OrnsteinUhlenbeckTradingLevelsExponential::optimalEntryLower(
    const double& d_star, const double& b_star, const double& stop_loss,
    const double& r, const double& c) const {
    // We need deep copies of the model and optimizer pointers to initialise the
    // params instance. This is because GSL requires a pointer to void and we
    // cannot use smart pointers with much benefit here. So we create deep
    // copies and then free that memory in the destructor of the
    // EntryLevelStopLossParams struct.
    void* params = new EntryLevelStopLossParams{
        newOptimizer(), newHittingTimeKernel(), b_star, stop_loss, r, c};
    ModelFunc fn = funcOptimalMeanReversionStopLossA;
    double value{0.0};
    try {
        double lower = stop_loss;
        double upper = d_star;
        value = brentSolver(fn, params, lower, upper);
    } catch (const std::exception& e) {
        std::cout
            << "Exception " << e.what()
            << " caught in "
               "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntryLower "
               "with stop loss."
            << std::endl;
        delete static_cast<EntryLevelStopLossParams*>(params);
        params = nullptr;
        throw;
    }
    // Then cast the void pointer back to the original type and free the memory.
    delete static_cast<EntryLevelStopLossParams*>(params);
    params = nullptr;
    return value;
}
const double OrnsteinUhlenbeckTradingLevelsExponential::optimalEntry(
    const double& b_star, const double& stop_loss, const double& r,
    const double& c) const {
    void* params = new EntryLevelStopLossParams{
        newOptimizer(), newHittingTimeKernel(), b_star, stop_loss, r, c};
    ModelFunc fn = funcOptimalMeanReversionStopLossD;
    double value{0.0};
    try {
        double lower = stop_loss;
        double upper = b_star;
        value = brentSolver(fn, params, lower, upper);
    } catch (const std::exception& e) {
        std::cout
            << "Exception " << e.what()
            << " caught in "
               "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntryLower "
               "with stop loss."
            << std::endl;
        delete static_cast<EntryLevelStopLossParams*>(params);
        params = nullptr;
        throw;
    }
    // Then cast the void pointer back to the original type and free the memory.
    delete static_cast<EntryLevelStopLossParams*>(params);
    params = nullptr;
    return value;
}
const double OrnsteinUhlenbeckTradingLevelsExponential::optimalEntry(
    const double& b_star, const double& r, const double& c) const {
    // We need deep copies of the model and optimizer pointers to initialise the
    // params instance. This is because GSL requires a pointer to void and we
    // cannot use smart pointers with much benefit here. So we create deep
    // copies and then free that memory in the destructor of the
    // EntryLevelParams struct.
    void* params = new EntryLevelParams{newOptimizer(), newHittingTimeKernel(),
                                        b_star, r, c};
    ModelFunc fn = funcOptimalMeanReversionD;
    double value{0.0};
    try {
        double lower = optimalEntryLowerBound();
        double upper = b_star;
        value = brentSolver(fn, params, lower, upper);
    } catch (const std::exception& e) {
        std::cout << "Exception " << e.what()
                  << " caught in "
                     "OrnsteinUhlenbeckTradingLevelsExponential::optimalEntry "
                     "without stop loss."
                  << std::endl;
        delete static_cast<EntryLevelParams*>(params);
        params = nullptr;
        throw;
    }
    // Then cast the void pointer back to the original type and free the memory.
    delete static_cast<EntryLevelParams*>(params);
    params = nullptr;
    return value;
}
