#include "stochastic_models/trading/optimal_mean_reversion.h"

#include <iostream>
#include <stdexcept>

#include "stochastic_models/numeric_utils/differentiation.h"
#include "stochastic_models/numeric_utils/integration.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"
#include "stochastic_models/trading/trading_levels_params.h"
OptimalMeanReversionParams::~OptimalMeanReversionParams() {
    delete hitting_time_kernel;
    hitting_time_kernel = nullptr;
}
double funcOptimalMeanReversionF(double x, void* params) {
    struct OptimalMeanReversionParams* p =
        static_cast<OptimalMeanReversionParams*>(params);
    double value = p->hitting_time_kernel->optimalTradingFCore(p->x, x, p->r);
    return value;
}
double valueFuncStopLoss(double x, void* params) {
    struct EntryLevelStopLossParams* p =
        static_cast<EntryLevelStopLossParams*>(params);
    double value = p->optimizer->V(p->hitting_time_kernel, x, p->b_star,
                                   p->stop_loss, p->r, p->c);
    return value;
}
double valueFunc(double x, void* params) {
    struct EntryLevelParams* p = static_cast<EntryLevelParams*>(params);
    double value =
        p->optimizer->V(p->hitting_time_kernel, x, p->b_star, p->r, p->c);
    return value;
}
double funcIntegrateF(double x, void* params) {
    struct ExitLevelParams* p = static_cast<ExitLevelParams*>(params);
    double value = p->optimizer->F(p->hitting_time_kernel, x, p->r, p->c);
    return value;
}
double funcIntegrateG(double x, void* params) {
    struct ExitLevelParams* p = static_cast<ExitLevelParams*>(params);
    double value = p->optimizer->G(p->hitting_time_kernel, x, p->r, p->c);
    return value;
}
double funcOptimalMeanReversionG(double x, void* params) {
    struct OptimalMeanReversionParams* p =
        static_cast<OptimalMeanReversionParams*>(params);
    double value = p->hitting_time_kernel->optimalTradingGCore(p->x, x, p->r);
    return value;
}
const OptimalMeanReversion* OptimalMeanReversion::clone() const {
    return new OptimalMeanReversion(*this);
}
const double OptimalMeanReversion::F(
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel, const double& x,
    const double& r, const double& c) const {
    // First create a deep copy of the model pointer and copy the contents of
    // model into the temporary location. This is because we are going to
    // free that memory in the destructor of the OptimalMeanReversionParams
    // struct.
    const HittingTimeOrnsteinUhlenbeck* temp_ptr = hitting_time_kernel->clone();

    // Create new instance of the OptimalMeanReversionParams struct as a void
    // pointer as this is required by the GSL numerical integration function.
    void* params = new OptimalMeanReversionParams{temp_ptr, x, r};
    ModelFunc fn = funcOptimalMeanReversionF;
    double lower = 0;
    double value = semiInfiniteIntegrationUpper(fn, params, lower);

    // Then cast the void pointer back to the original type and free the memory.
    OptimalMeanReversionParams* ptr =
        static_cast<OptimalMeanReversionParams*>(params);
    delete ptr;
    ptr = nullptr;
    params = nullptr;
    temp_ptr = nullptr;

    return value;
}
const double OptimalMeanReversion::G(
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel, const double& x,
    const double& r, const double& c) const {
    // First create a deep copy of the model pointer and copy the contents of
    // model into the temporary location. This is because we are going to free
    // that memory in the destructor of the OptimalMeanReversionParams struct.
    const HittingTimeOrnsteinUhlenbeck* temp_ptr = hitting_time_kernel->clone();

    // Create new instance of the OptimalMeanReversionParams struct as a void
    // pointer as this is required by the GSL numerical integration function.
    void* params = new OptimalMeanReversionParams{temp_ptr, x, r};
    ModelFunc fn = funcOptimalMeanReversionG;
    double value = 0.0;
    try {
        double lower = 0;
        value = semiInfiniteIntegrationUpper(fn, params, lower);
    } catch (const std::exception& e) {
        std::cout << "Exception " << e.what()
                  << " caught in OptimalMeanReversion::G." << std::endl;
        struct OptimalMeanReversionParams* p =
            static_cast<OptimalMeanReversionParams*>(params);
        delete p;
        p = nullptr;
        params = nullptr;
        temp_ptr = nullptr;
        throw;
    }
    // Then cast the void pointer back to the original type and free the memory.
    OptimalMeanReversionParams* ptr =
        static_cast<OptimalMeanReversionParams*>(params);
    delete ptr;
    ptr = nullptr;
    params = nullptr;
    temp_ptr = nullptr;

    return value;
}
const double OptimalMeanReversion::L_star(
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel, const double& r,
    const double& c) const {
    return hitting_time_kernel->optimalTradingLCore(r, c);
}
const double OptimalMeanReversion::b(
    const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel, const double& r,
    const double& c) const {
    ModelFunc fn = funcIntegrateF;
    const double result =
        OptimalMeanReversion::F(hitting_time_kernel, value, r, c) -
        (value - c) *
            instantaneousDifferential(fn, hitting_time_kernel, value, r, c);
    return result;
}
const double OptimalMeanReversion::b(
    const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& stop_loss, const double& r, const double& c) const {
    ModelFunc fnF = funcIntegrateF;
    ModelFunc fnG = funcIntegrateG;

    const double bMinusC = value - c;
    const double lMinusC = stop_loss - c;
    const double gL =
        OptimalMeanReversion::G(hitting_time_kernel, stop_loss, r, c);
    const double gB = OptimalMeanReversion::G(hitting_time_kernel, value, r, c);
    const double fL =
        OptimalMeanReversion::F(hitting_time_kernel, stop_loss, r, c);
    const double fB = OptimalMeanReversion::F(hitting_time_kernel, value, r, c);

    const double result =
        (((lMinusC * gB) - (bMinusC * gL)) *
         instantaneousDifferential(fnF, hitting_time_kernel, value, r, c)) +
        (((bMinusC * fL) - (lMinusC * fB)) *
         instantaneousDifferential(fnG, hitting_time_kernel, value, r, c)) -
        ((gB * fL) - (gL * fB));
    return result;
}
const double OptimalMeanReversion::d(
    const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star, const double& stop_loss, const double& r,
    const double& c) const {
    ModelFunc funcG = funcIntegrateG;
    ModelFunc funcV = valueFuncStopLoss;
    const double vPrimeD = instantaneousDifferential(
        funcV, hitting_time_kernel, value, b_star, stop_loss, r, c);
    const double result =
        (OptimalMeanReversion::G(hitting_time_kernel, value, r, c) *
         (vPrimeD - 1)) -
        (instantaneousDifferential(funcG, hitting_time_kernel, value, b_star, r,
                                   c) *
         (OptimalMeanReversion::V(hitting_time_kernel, value, b_star, stop_loss,
                                  r, c) -
          value - c));

    return result;
}
const double OptimalMeanReversion::d(
    const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star, const double& r, const double& c) const {
    ModelFunc funcG = funcIntegrateG;
    ModelFunc funcV = valueFunc;
    const double result =
        (OptimalMeanReversion::G(hitting_time_kernel, value, r, c) *
         (instantaneousDifferential(funcV, hitting_time_kernel, value, b_star,
                                    r, c) -
          1)) -
        (instantaneousDifferential(funcG, hitting_time_kernel, value, b_star, r,
                                   c) *
         (OptimalMeanReversion::V(hitting_time_kernel, value, b_star, r, c) -
          value - c));

    return result;
}
const double OptimalMeanReversion::a(
    const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star, const double& r, const double& c) const {
    // There is no mathematical definition for the optimal exit level of
    // a mean reversion model with no stop loss parameter.
    throw std::logic_error(
        "OptimalMeanReversion does not implement function "
        "a with no stop_loss parameter");
}
const double OptimalMeanReversion::a(
    const double& value,
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& b_star, const double& stop_loss, const double& r,
    const double& c) const {
    ModelFunc funcF = funcIntegrateF;
    ModelFunc funcV = valueFuncStopLoss;
    const double vPrimeA = instantaneousDifferential(
        funcV, hitting_time_kernel, value, b_star, stop_loss, r, c);
    const double result =
        (OptimalMeanReversion::F(hitting_time_kernel, value, r, c) *
         (vPrimeA - 1)) -
        (instantaneousDifferential(funcF, hitting_time_kernel, value, b_star, r,
                                   c) *
         (OptimalMeanReversion::V(hitting_time_kernel, value, b_star, stop_loss,
                                  r, c) -
          value - c));

    return result;
}
const double OptimalMeanReversion::V(
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel, const double& x,
    const double& b_star, const double& r, const double& c) const {
    if (x < b_star) {
        return (b_star - c) * F(hitting_time_kernel, x, r, c) /
               F(hitting_time_kernel, b_star, r, c);
    } else {
        return x - c;
    }
}
const double OptimalMeanReversion::V(
    const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel, const double& x,
    const double& b_star, const double& stop_loss, const double& r,
    const double& c) const {
    if ((b_star > x) && (x > stop_loss)) {
        ModelFunc fnF = funcIntegrateF;
        ModelFunc fnG = funcIntegrateG;

        const double bMinusC = b_star - c;
        const double lMinusC = stop_loss - c;
        const double gL =
            OptimalMeanReversion::G(hitting_time_kernel, stop_loss, r, c);
        const double gBstar =
            OptimalMeanReversion::G(hitting_time_kernel, b_star, r, c);
        const double fL =
            OptimalMeanReversion::F(hitting_time_kernel, stop_loss, r, c);
        const double fBstar =
            OptimalMeanReversion::F(hitting_time_kernel, b_star, r, c);

        const double C = ((bMinusC * gL) - (lMinusC * gBstar)) /
                         ((fBstar * gL) - (fL * gBstar));
        const double D = ((lMinusC * fBstar) - (bMinusC * fL)) /
                         ((fBstar * gL) - (fL * gBstar));

        return C * F(hitting_time_kernel, x, r, c) +
               D * G(hitting_time_kernel, x, r, c);
    } else {
        return x - c;
    }
}
const double OptimalMeanReversion::instantaneousDifferential(
    ModelFunc fn, const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x, const double& b_star, const double& stop_loss,
    const double& r, const double& c) const {
    // Create a deep copy of the model and optimizer pointers to ensure that
    // cleanup does not delete heap memory that is not owned by this scope.
    const HittingTimeOrnsteinUhlenbeck* temp_kernel =
        hitting_time_kernel->clone();
    const OptimalMeanReversion* temp_optimizer = clone();

    // Create copy to not modify the original value and
    // not provide to function that takes a non-const.
    double x_copied = x;

    // Adaptive differentiation function.
    void* params = new EntryLevelStopLossParams{
        temp_optimizer, temp_kernel, b_star, stop_loss, r, c};

    try {
        double value = adaptiveCentralDifferentiation(fn, params, x_copied);
        struct EntryLevelStopLossParams* p =
            static_cast<EntryLevelStopLossParams*>(params);
        delete p;
        p = nullptr;
        params = nullptr;
        return value;
    } catch (std::exception& e) {
        std::cout
            << "Exception " << e.what()
            << " caught in OptimalMeanReversion::instantaneousDifferential."
            << std::endl;
        struct EntryLevelStopLossParams* p =
            static_cast<EntryLevelStopLossParams*>(params);
        delete p;
        p = nullptr;
        params = nullptr;
        throw;
    }
}
const double OptimalMeanReversion::instantaneousDifferential(
    ModelFunc fn, const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x, const double& b_star, const double& r,
    const double& c) const {
    // Create a deep copy of the model and optimizer pointers to ensure that
    // cleanup does not delete heap memory that is not owned by this scope.
    const HittingTimeOrnsteinUhlenbeck* temp_kernel =
        hitting_time_kernel->clone();
    const OptimalMeanReversion* temp_optimizer = clone();

    // Create copy to not modify the original value and
    // not provide to function that takes a non-const.
    double x_copied = x;

    // Adaptive differentiation function.
    void* params =
        new EntryLevelParams{temp_optimizer, temp_kernel, b_star, r, c};

    try {
        double value = adaptiveCentralDifferentiation(fn, params, x_copied);
        struct EntryLevelParams* p = static_cast<EntryLevelParams*>(params);
        delete p;
        p = nullptr;
        params = nullptr;
        return value;
    } catch (std::exception& e) {
        std::cout
            << "Exception " << e.what()
            << " caught in OptimalMeanReversion::instantaneousDifferential."
            << std::endl;
        struct EntryLevelParams* p = static_cast<EntryLevelParams*>(params);
        delete p;
        p = nullptr;
        params = nullptr;
        throw;
    }
}
const double OptimalMeanReversion::instantaneousDifferential(
    ModelFunc fn, const HittingTimeOrnsteinUhlenbeck* hitting_time_kernel,
    const double& x, const double& r, const double& c) const {
    // Create a deep copy of the model and optimizer pointers to ensure that
    // cleanup does not delete heap memory that is not owned by this scope.
    const HittingTimeOrnsteinUhlenbeck* temp_kernel =
        hitting_time_kernel->clone();
    const OptimalMeanReversion* temp_optimizer = clone();

    // Create copy to not modify the original value and
    // not provide to function that takes a non-const.
    double x_copied = x;

    // Adaptive differentiation function.
    void* params = new ExitLevelParams{temp_optimizer, temp_kernel, r, c};

    try {
        double value = adaptiveCentralDifferentiation(fn, params, x_copied);
        struct ExitLevelParams* p = static_cast<ExitLevelParams*>(params);
        delete p;
        p = nullptr;
        params = nullptr;
        return value;
    } catch (std::exception& e) {
        std::cout
            << "Exception " << e.what()
            << " caught in OptimalMeanReversion::instantaneousDifferential."
            << std::endl;
        struct ExitLevelParams* p = static_cast<ExitLevelParams*>(params);
        delete p;
        p = nullptr;
        params = nullptr;
        throw;
    }
}
