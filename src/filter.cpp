#include "stochastic_models/kalman_filter/filter.h"

#include <utility>

#include "stochastic_models/kalman_filter/states_exceptions.h"

KineticComponents::KineticComponents(const FilterSystemDimensions &dimensions)
    : filter_state(dimensions) {}
void KineticComponents::setFilterState(const KcaStates &state) {
    filter_state = state;
}
const KcaStates KineticComponents::getFilterState() const {
    return filter_state;
}
const bool &KineticComponents::isInitialised() const {
    return filter_state.isInitialised();
}
const bool &KineticComponents::isPriorStateValid() const {
    return filter_state.arePriorsValid();
}
const std::vector<double> KineticComponents::getCurrentState() const {
    return filter_state.getCurrentStateMeanVector();
}
void KineticComponents::initialiseFilter(const std::vector<double> &data_series,
                                         const double &h, const double &q) {
    filter_state.setInitialState(data_series, h, q);
}
void KineticComponents::updatePriors() {
    try {
        filter_state.updatePredictedState();
    } catch (const filter_uninitialised &exc) {
        std::string message =
            "Unhandled error when updating the prior kinetic components "
            "state: " +
            std::string(exc.what());
        throw filter_uninitialised(message);
    }
}
void KineticComponents::updatePosteriors(const double &observation,
                                         const double &innovation_sigma) {
    try {
        filter_state.updateCurrentState(observation, innovation_sigma);
    } catch (const filter_invalid_operation &exc) {
        std::string message =
            "Unhandled error when updating the posterior kinetic components "
            "state: " +
            std::string(exc.what());
        throw filter_invalid_operation("");
    }
}
