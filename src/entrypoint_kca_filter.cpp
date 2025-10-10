/**
 * @file entrypoint_kca_filter.cpp
 * @brief Entry point module for Kinetic Components Analysis initialization and
 * update steps.
 *
 */
#include <string>
#include <vector>

#include "stochastic_models/entrypoints/kca_filter.h"
#include "stochastic_models/kalman_filter/adapters.h"
#include "stochastic_models/kalman_filter/kca.h"

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

    // Create JSON adapter to handle serialisation and deserialisation of
    // the internal state provided.
    KcaStatesJsonAdapter adapter;

    // Set the internal state of the filter with the state matrices /
    // vectors provided.
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
