#include "stochastic_models/trading/trading_levels.h"

#include <algorithm>
#include <iostream>

#include "stochastic_models/numeric_utils/helpers.h"
#include "stochastic_models/numeric_utils/solvers.h"
ExitLevelStopLossParams::~ExitLevelStopLossParams() {
    delete optimizer;
    optimizer = nullptr;
    delete hitting_time_kernel;
    hitting_time_kernel = nullptr;
}
ExitLevelParams::~ExitLevelParams() {
    delete optimizer;
    optimizer = nullptr;
    delete hitting_time_kernel;
    hitting_time_kernel = nullptr;
}
EntryLevelStopLossParams::~EntryLevelStopLossParams() {
    delete optimizer;
    optimizer = nullptr;
    delete hitting_time_kernel;
    hitting_time_kernel = nullptr;
}
EntryLevelParams::~EntryLevelParams() {
    delete optimizer;
    optimizer = nullptr;
    delete hitting_time_kernel;
    hitting_time_kernel = nullptr;
}
double funcOptimalMeanReversionA(double x, void* params) {
    struct EntryLevelParams* p = static_cast<EntryLevelParams*>(params);
    return p->optimizer->a(x, p->hitting_time_kernel, p->b_star, p->r, p->c);
}
double funcOptimalMeanReversionStopLossA(double x, void* params) {
    struct EntryLevelStopLossParams* p =
        static_cast<EntryLevelStopLossParams*>(params);
    return p->optimizer->a(x, p->hitting_time_kernel, p->b_star, p->stop_loss,
                           p->r, p->c);
}
double funcOptimalMeanReversionStopLossD(double x, void* params) {
    struct EntryLevelStopLossParams* p =
        static_cast<EntryLevelStopLossParams*>(params);
    return p->optimizer->d(x, p->hitting_time_kernel, p->b_star, p->stop_loss,
                           p->r, p->c);
}
double funcOptimalMeanReversionStopLossB(double x, void* params) {
    struct ExitLevelStopLossParams* p =
        static_cast<ExitLevelStopLossParams*>(params);
    return p->optimizer->b(x, p->hitting_time_kernel, p->stop_loss, p->r, p->c);
}
double funcOptimalMeanReversionD(double x, void* params) {
    struct EntryLevelParams* p = static_cast<EntryLevelParams*>(params);
    return p->optimizer->d(x, p->hitting_time_kernel, p->b_star, p->r, p->c);
}
double funcOptimalMeanReversionB(double x, void* params) {
    struct ExitLevelParams* p = static_cast<ExitLevelParams*>(params);
    return p->optimizer->b(x, p->hitting_time_kernel, p->r, p->c);
}
const double TradingLevels::optimalExit(
    const OptimalTrading* optimizer, const StochasticModel* model,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& stop_loss, const double& r, const double& c) const {
    // Create core model and optimal mean reversion instances.
    // The lower bound is the maximum of the cost of trading and the L function.
    double lower = std::max(optimizer->L_star(hitting_time_kernel, r, c), c);
    double upper = upperSolverBound(model);

    // First create a deep copy of the model pointer and copy the contents of
    // model into the temporary location. This is because we are going to
    // free that memory in the destructor of the ExitLevelStopLossParams
    // struct.
    const HittingTimeOrnsteinUhlenbeck* temp_kernel =
        hitting_time_kernel->clone();
    const OptimalTrading* temp_optimizer = optimizer->clone();

    void* params = new ExitLevelStopLossParams{temp_optimizer, temp_kernel,
                                               stop_loss, r, c};
    ModelFunc fn = funcOptimalMeanReversionStopLossB;
    double value = brentSolver(fn, params, lower, upper);

    // Then cast the void pointer back to the original type and free the memory.
    ExitLevelStopLossParams* ptr =
        static_cast<ExitLevelStopLossParams*>(params);
    delete ptr;
    ptr = nullptr;
    params = nullptr;
    temp_optimizer = nullptr;
    temp_kernel = nullptr;
    return value;
}
const double TradingLevels::optimalExit(
    const OptimalTrading* optimizer, const StochasticModel* model,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel, const double& r,
    const double& c) const {
    // Create core model and optimal mean reversion instances.
    // The lower bound is the maximum of the cost of trading and the L function.
    double lower = std::max(optimizer->L_star(hitting_time_kernel, r, c), c);
    double upper = upperSolverBound(model);

    // First create a deep copy of the model pointer and copy the contents of
    // model into the temporary location. This is because we are going to
    // free that memory in the destructor of the ExitLevelParams
    // struct.
    const HittingTimeOrnsteinUhlenbeck* temp_kernel =
        hitting_time_kernel->clone();
    const OptimalTrading* temp_optimizer = optimizer->clone();

    void* params = new ExitLevelParams{temp_optimizer, temp_kernel, r, c};
    ModelFunc fn = funcOptimalMeanReversionB;

    double value{0.0};
    try {
        value = brentSolver(fn, params, lower, upper);
    } catch (const std::runtime_error& e) {
        std::cout << "Exception " << e.what()
                  << " caught in TradingLevels::optimalExit without stop loss."
                  << std::endl;
        delete static_cast<ExitLevelParams*>(params);
        params = nullptr;
        temp_optimizer = nullptr;
        temp_kernel = nullptr;
        throw;
    }
    // Then cast the void pointer back to the original type and free the memory.
    ExitLevelParams* ptr = static_cast<ExitLevelParams*>(params);
    delete ptr;
    ptr = nullptr;
    params = nullptr;
    temp_optimizer = nullptr;
    temp_kernel = nullptr;
    return value;
}
const double TradingLevels::optimalEntryLower(
    const OptimalTrading* optimizer, const StochasticModel* model,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& d_star, const double& b_star, const double& r,
    const double& c) const {
    // Create core model and optimal mean reversion instances.
    double lower = lowerSolverBound(model);
    double upper = d_star;

    // First create a deep copy of the model pointer and copy the contents
    // of model into the temporary location. This is because we are going to
    // free that memory in the destructor of the EntryLevelParams
    // struct.
    const HittingTimeOrnsteinUhlenbeck* temp_kernel =
        hitting_time_kernel->clone();
    const OptimalTrading* temp_optimizer = optimizer->clone();

    void* params =
        new EntryLevelParams{temp_optimizer, temp_kernel, b_star, r, c};
    ModelFunc fn = funcOptimalMeanReversionA;
    double value = brentSolver(fn, params, lower, upper);

    // Then cast the void pointer back to the original type and free the memory.
    EntryLevelParams* ptr = static_cast<EntryLevelParams*>(params);
    delete ptr;
    ptr = nullptr;
    params = nullptr;
    temp_optimizer = nullptr;
    temp_kernel = nullptr;
    return value;
}
const double TradingLevels::optimalEntryLower(
    const OptimalTrading* optimizer, const StochasticModel* model,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& d_star, const double& b_star, const double& stop_loss,
    const double& r, const double& c) const {
    // Create core model and optimal mean reversion instances.
    double lower = stop_loss;
    double upper = d_star;

    // First create a deep copy of the model pointer and copy the contents
    // of model into the temporary location. This is because we are going to
    // free that memory in the destructor of the EntryLevelParams
    // struct.
    const HittingTimeOrnsteinUhlenbeck* temp_kernel =
        hitting_time_kernel->clone();
    const OptimalTrading* temp_optimizer = optimizer->clone();

    void* params = new EntryLevelStopLossParams{
        temp_optimizer, temp_kernel, b_star, stop_loss, r, c};
    ModelFunc fn = funcOptimalMeanReversionStopLossA;
    double value = brentSolver(fn, params, lower, upper);

    // Then cast the void pointer back to the original type and free the memory.
    EntryLevelStopLossParams* ptr =
        static_cast<EntryLevelStopLossParams*>(params);
    delete ptr;
    ptr = nullptr;
    params = nullptr;
    temp_optimizer = nullptr;
    temp_kernel = nullptr;
    return value;
}
const double TradingLevels::optimalEntry(
    const OptimalTrading* optimizer, const StochasticModel* model,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star, const double& stop_loss, const double& r,
    const double& c) const {
    // Create core model and optimal mean reversion instances.
    double lower = stop_loss;
    double upper = b_star;

    // First create a deep copy of the model pointer and copy the contents
    // of model into the temporary location. This is because we are going to
    // free that memory in the destructor of the EntryLevelParams
    // struct.
    const HittingTimeOrnsteinUhlenbeck* temp_kernel =
        hitting_time_kernel->clone();
    const OptimalTrading* temp_optimizer = optimizer->clone();

    void* params = new EntryLevelStopLossParams{
        temp_optimizer, temp_kernel, b_star, stop_loss, r, c};
    ModelFunc fn = funcOptimalMeanReversionStopLossD;
    double value = brentSolver(fn, params, lower, upper);

    // Then cast the void pointer back to the original type and free the memory.
    EntryLevelStopLossParams* ptr =
        static_cast<EntryLevelStopLossParams*>(params);
    delete ptr;
    ptr = nullptr;
    params = nullptr;
    temp_optimizer = nullptr;
    temp_kernel = nullptr;
    return value;
}
const double TradingLevels::optimalEntry(
    const OptimalTrading* optimizer, const StochasticModel* model,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star, const double& r, const double& c) const {
    // Create core model and optimal mean reversion instances.
    double lower = lowerSolverBound(model);
    double upper = b_star;

    // First create a deep copy of the model pointer and copy the contents
    // of model into the temporary location. This is because we are going to
    // free that memory in the destructor of the EntryLevelParams
    // struct.
    const HittingTimeOrnsteinUhlenbeck* temp_kernel =
        hitting_time_kernel->clone();
    const OptimalTrading* temp_optimizer = optimizer->clone();

    void* params =
        new EntryLevelParams{temp_optimizer, temp_kernel, b_star, r, c};
    ModelFunc fn = funcOptimalMeanReversionD;
    double value = brentSolver(fn, params, lower, upper);

    // Then cast the void pointer back to the original type and free the memory.
    EntryLevelParams* ptr = static_cast<EntryLevelParams*>(params);
    delete ptr;
    ptr = nullptr;
    params = nullptr;
    temp_optimizer = nullptr;
    temp_kernel = nullptr;
    return value;
}
