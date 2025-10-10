#ifndef KCA_FILTER_H
#define KCA_FILTER_H
#include <string>
#include <vector>

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

#endif  // KCA_FILTER_H
