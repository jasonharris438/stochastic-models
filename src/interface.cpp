/**
 * @file interface.cpp
 * @brief Interface module to the stochastic models library.
 * @version 0.1.1
 *
 */
#include "stochastic_models/hitting_times/hitting_time_density.h"
#include "stochastic_models/interface/stochastic_models.h"
#include "stochastic_models/kalman_filter/adapters.h"
#include "stochastic_models/likelihood/ornstein_uhlenbeck_likelihood.h"
#include "stochastic_models/sde/ornstein_uhlenbeck.h"
#include "stochastic_models/trading/exponential_mean_reversion.h"
#include "stochastic_models/trading/optimal_mean_reversion.h"
#include "stochastic_models/trading/trading_levels.h"

const double optimalExitLevel(const double mu, const double alpha,
                              const double sigma, const double stop_loss,
                              const double r, const double c) {
    // Create core model, optimal mean reversion, and trading levels instances.
    // Model containing SDE equation and parameters.
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(mu, alpha, sigma);

    // Object with trading levels equation components interface.
    OptimalMeanReversion *optimizer = new OptimalMeanReversion();

    // Object with outer/controller layer logic for trading levels.
    TradingLevels tradingLevels;

    // Object with trading levels kernel layer interface.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);

    // Calculate b*.
    const double value = tradingLevels.optimalExit(
        optimizer, model, hitting_time_kernel, stop_loss, r, c);
    delete model;
    delete hitting_time_kernel;
    delete optimizer;
    return value;
}
const double optimalExitLevel(const double mu, const double alpha,
                              const double sigma, const double r,
                              const double c) {
    // Create core model, optimal mean reversion, and trading levels instances.
    // Model containing SDE equation and parameters.
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(mu, alpha, sigma);

    // Object with trading levels equation components interface.
    OptimalMeanReversion *optimizer = new OptimalMeanReversion();

    // Object with outer/controller layer logic for trading levels.
    TradingLevels tradingLevels;

    // Object with trading levels kernel layer interface.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);

    // Calculate b*.
    const double value =
        tradingLevels.optimalExit(optimizer, model, hitting_time_kernel, r, c);
    delete model;
    delete hitting_time_kernel;
    delete optimizer;
    return value;
}
const double optimalExitLevelExponential(const double mu, const double alpha,
                                         const double sigma, const double r,
                                         const double c) {
    // Create core model, optimal mean reversion, and trading levels instances.
    // Model containing SDE equation and parameters.
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(mu, alpha, sigma);

    // Object with trading levels equation components interface.
    ExponentialMeanReversion *optimizer = new ExponentialMeanReversion();

    // Object with outer/controller layer logic for trading levels.
    TradingLevels tradingLevels;

    // Object with trading levels kernel layer interface.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);

    double value{0.0};
    try {
        // Calculate b*.
        value = tradingLevels.optimalExit(optimizer, model, hitting_time_kernel,
                                          r, c);
    } catch (const std::runtime_error &e) {
        std::cout << "Exception " << e.what()
                  << " caught in optimalExitLevelExponential." << std::endl;
        delete model;
        delete hitting_time_kernel;
        delete optimizer;
        throw;
    }
    delete model;
    delete hitting_time_kernel;
    delete optimizer;
    return value;
}
const double optimalEntryLevel(const double b_star, const double mu,
                               const double alpha, const double sigma,
                               const double r, const double c) {
    // Create core model, optimal mean reversion, and trading levels instances.
    // Model containing SDE equation and parameters.
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(mu, alpha, sigma);

    // Object with trading levels equation components interface.
    OptimalMeanReversion *optimizer = new OptimalMeanReversion();

    // Object with outer/controller layer logic for trading levels.
    TradingLevels tradingLevels;

    // Object with trading levels kernel layer interface.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);

    // Calculate d*.
    const double value = tradingLevels.optimalEntry(
        optimizer, model, hitting_time_kernel, b_star, r, c);
    delete model;
    delete hitting_time_kernel;
    delete optimizer;
    return value;
}
const double optimalEntryLevelExponential(const double b_star, const double mu,
                                          const double alpha,
                                          const double sigma, const double r,
                                          const double c) {
    // Create core model, optimal mean reversion, and trading levels instances.
    // Model containing SDE equation and parameters.
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(mu, alpha, sigma);

    // Object with trading levels equation components interface.
    ExponentialMeanReversion *optimizer = new ExponentialMeanReversion();

    // Object with outer/controller layer logic for trading levels.
    TradingLevels tradingLevels;

    // Object with trading levels kernel layer interface.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);

    // Calculate d*.
    const double value = tradingLevels.optimalEntry(
        optimizer, model, hitting_time_kernel, b_star, r, c);
    delete model;
    delete hitting_time_kernel;
    delete optimizer;
    return value;
}
const double optimalEntryLevelLowerExponential(
    const double d_star, const double b_star, const double mu,
    const double alpha, const double sigma, const double r, const double c) {
    // Create core model, optimal mean reversion, and trading levels instances.
    // Model containing SDE equation and parameters.
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(mu, alpha, sigma);

    // Object with trading levels equation components interface.
    ExponentialMeanReversion *optimizer = new ExponentialMeanReversion();

    // Object with outer/controller layer logic for trading levels.
    TradingLevels tradingLevels;

    // Object with trading levels kernel layer interface.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);

    // Calculate d*.
    const double value = tradingLevels.optimalEntryLower(
        optimizer, model, hitting_time_kernel, d_star, b_star, r, c);
    delete model;
    delete hitting_time_kernel;
    delete optimizer;
    return value;
}
const double optimalEntryLevelLower(const double d_star, const double b_star,
                                    const double mu, const double alpha,
                                    const double sigma, const double stop_loss,
                                    const double r, const double c) {
    // Create core model, optimal mean reversion, and trading levels instances.
    // Model containing SDE equation and parameters.
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(mu, alpha, sigma);

    // Object with trading levels equation components interface.
    OptimalMeanReversion *optimizer = new OptimalMeanReversion();

    // Object with outer/controller layer logic for trading levels.
    TradingLevels tradingLevels;

    // Object with trading levels kernel layer interface.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);

    // Calculate d*.
    const double value = tradingLevels.optimalEntryLower(
        optimizer, model, hitting_time_kernel, d_star, b_star, stop_loss, r, c);
    delete model;
    delete hitting_time_kernel;
    delete optimizer;
    return value;
}
const double optimalEntryLevel(const double b_star, const double mu,
                               const double alpha, const double sigma,
                               const double stop_loss, const double r,
                               const double c) {
    // Create core model, optimal mean reversion, and trading levels instances.
    // Model containing SDE equation and parameters.
    OrnsteinUhlenbeckModel *model =
        new OrnsteinUhlenbeckModel(mu, alpha, sigma);

    // Object with trading levels equation components interface.
    OptimalMeanReversion *optimizer = new OptimalMeanReversion();

    // Object with outer/controller layer logic for trading levels.
    TradingLevels tradingLevels;

    // Object with trading levels kernel layer interface.
    HittingTimeOrnsteinUhlenbeck *hitting_time_kernel =
        new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);

    // Calculate d*.
    const double value = tradingLevels.optimalEntry(
        optimizer, model, hitting_time_kernel, b_star, stop_loss, r, c);
    delete model;
    delete hitting_time_kernel;
    delete optimizer;
    return value;
}
const double hittingTimeDensityOrnsteinUhlenbeck(double x, const double mu,
                                                 const double alpha,
                                                 const double sigma,
                                                 double first, double second) {
    // Create core model as pointer to void because that type is required by
    // GSL.
    void *model = new HittingTimeOrnsteinUhlenbeck(mu, alpha, sigma);

    // The ModelFunc is provided to determine which function to provide to GSL
    // for integration.
    ModelFunc fn = &integrateHittingTimeDensity;

    // Calculate the hitting time density at point x given the arguments first
    // and second.
    const double value = hittingTimeDensity(x, fn, model, first, second);

    // Then need to cast model back to OrnsteinUhlenbeckModel pointer to clean
    // up.
    HittingTimeOrnsteinUhlenbeck *ptr =
        static_cast<HittingTimeOrnsteinUhlenbeck *>(model);
    delete ptr;
    ptr = nullptr;
    model = nullptr;
    return value;
}
const std::unordered_map<std::string, const double>
ornsteinUhlenbeckMaximumLikelihood(const std::vector<double> vec) {
    // Generate likelihood calculator and generate estimate of mu, alpha, and
    // sigma.
    OrnsteinUhlenbeckLikelihood *likelihood = new OrnsteinUhlenbeckLikelihood();
    const std::unordered_map<std::string, const double> params{
        likelihood->calculate(vec)};
    delete likelihood;
    return params;
}
const std::string getInitializedKcaState(const std::vector<double> data_series,
                                         const double h, const double q,
                                         const std::string system_dimensions) {
    // Deserialize the initial dimensions of the filter system that were
    // provided in the initial state.
    const FilterSystemDimensionsJsonAdapter dimensions_adapter;
    const FilterSystemDimensions dimensions =
        dimensions_adapter.deserialize(system_dimensions);

    // Create kinetic components object.
    KineticComponents kinetic_components = KineticComponents{dimensions};

    // Initialize the filter with the provided data series and parameters.
    kinetic_components.initialiseFilter(data_series, h, q);

    // Return the internal state of the filter.
    const KcaStates internal_state = kinetic_components.getFilterState();

    const KcaStatesJsonAdapter state_adapter;
    return state_adapter.serialize(internal_state);
}
const std::string getUpdatedKcaState(const std::string state,
                                     const std::string system_dimensions,
                                     const double observation,
                                     const double innovation_sigma) {
    // Deserialize the initial dimensions of the filter system that were
    // provided in the initial state.
    const FilterSystemDimensionsJsonAdapter dimensions_adapter;
    const FilterSystemDimensions dimensions =
        dimensions_adapter.deserialize(system_dimensions);

    // Create JSON adapter to handle serialisation and deserialisation of the
    // internal state provided.
    KcaStatesJsonAdapter adapter;

    // Set the internal state of the filter with the state matrices / vectors
    // provided.
    KcaStates internal_state = adapter.deserialize(state, dimensions);

    // Create kinetic components object and set the filter state with the
    // provided internal state.
    KineticComponents kinetic_components = KineticComponents{dimensions};
    kinetic_components.setFilterState(internal_state);

    // Update priors and posteriors with the provided observation and
    // innovation sigma.
    kinetic_components.updatePriors();
    kinetic_components.updatePosteriors(observation, innovation_sigma);

    // Return the internal state of the filter.
    const KcaStates updated_state = kinetic_components.getFilterState();

    // serialize the internal state to a JSON string before returning.
    return adapter.serialize(updated_state);
}
