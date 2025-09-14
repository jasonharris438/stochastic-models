#include "stochastic_models/kalman_filter/adapters.h"

#include "stochastic_models/kalman_filter/states.h"
#include "stochastic_models/kalman_filter/states_exceptions.h"

// nlohmann suggested to improve the depth of error reporting from json objects.
// When enabled, exception messages contain a JSON Pointer to the JSON value
// that triggered the exception. This carries additional runtime overhead.
// https://json.nlohmann.me/home/exceptions/#extended-diagnostic-messages
#ifndef JSON_DIAGNOSTICS
#define JSON_DIAGNOSTICS 0
#endif
#include <nlohmann/json.hpp>
const std::vector<std::vector<double>>
KcaStatesJsonAdapter::copyBoostMatrixToVector(
    const boost::numeric::ublas::matrix<double> &boost_matrix) const {
    std::vector<std::vector<double>> result(
        boost_matrix.size1(), std::vector<double>(boost_matrix.size2()));

    // Copy the values from the matrix to the std::vector by row.
    for (int i{0}; i < boost_matrix.size1(); i++) {
        std::copy(row(boost_matrix, i).begin(), row(boost_matrix, i).end(),
                  result.at(i).begin());
    }
    return result;
}
const std::vector<double> KcaStatesJsonAdapter::copyBoostVectorToVector(
    const boost::numeric::ublas::vector<double> &boost_vector) const {
    std::vector<double> std_vector(boost_vector.size());
    std::copy(boost_vector.begin(), boost_vector.end(), std_vector.begin());
    return std_vector;
}
const FilterSystemDimensions FilterSystemDimensionsJsonAdapter::deserialize(
    const std::string &state) const {
    nlohmann::json json_obj = nlohmann::json::parse(state);

    try {
        FilterSystemDimensions dimensions;
        json_obj.at("state_mean_dimension")
            .get_to(dimensions.state_mean_dimension);
        json_obj.at("state_covariance_rows")
            .get_to(dimensions.state_covariance_rows);
        json_obj.at("state_covariance_columns")
            .get_to(dimensions.state_covariance_columns);
        json_obj.at("observation_matrix_rows")
            .get_to(dimensions.observation_matrix_rows);
        json_obj.at("observation_matrix_columns")
            .get_to(dimensions.observation_matrix_columns);
        json_obj.at("observation_covariance_rows")
            .get_to(dimensions.observation_covariance_rows);
        json_obj.at("observation_covariance_columns")
            .get_to(dimensions.observation_covariance_columns);
        json_obj.at("observation_offset").get_to(dimensions.observation_offset);

        return dimensions;
    } catch (const nlohmann::json::exception &exc) {
        throw json_parse_error(exc.what());
    }
}
const std::string FilterSystemDimensionsJsonAdapter::serialize(
    const FilterSystemDimensions &dimensions) const {
    nlohmann::json json_obj = {
        {"state_mean_dimension", dimensions.state_mean_dimension},
        {"state_covariance_rows", dimensions.state_covariance_rows},
        {"state_covariance_columns", dimensions.state_covariance_columns},
        {"observation_matrix_rows", dimensions.observation_matrix_rows},
        {"observation_matrix_columns", dimensions.observation_matrix_columns},
        {"observation_covariance_rows", dimensions.observation_covariance_rows},
        {"observation_covariance_columns",
         dimensions.observation_covariance_columns},
        {"observation_offset", dimensions.observation_offset}};
    return json_obj.dump();
}
const std::string KcaStatesJsonAdapter::serialize(
    const KcaStates &kca_states) const {
    nlohmann::json json_obj;
    json_obj["transition_matrix"] =
        copyBoostMatrixToVector(kca_states.getTransitionMatrix());
    json_obj["transition_covariance"] =
        copyBoostMatrixToVector(kca_states.getTransitionCovariance());
    json_obj["current_state_covariance"] =
        copyBoostMatrixToVector(kca_states.getCurrentStateCovariance());
    json_obj["current_state_mean"] =
        copyBoostVectorToVector(kca_states.getCurrentStateMean());
    json_obj["observation_matrix"] =
        copyBoostMatrixToVector(kca_states.getObservationMatrix());
    json_obj["observation_offset"] = kca_states.getObservationOffset();
    return json_obj.dump();
}
const KcaStates KcaStatesJsonAdapter::deserialize(
    const std::string &state, const FilterSystemDimensions &dimensions) const {
    nlohmann::json json_obj = nlohmann::json::parse(state);

    KcaStates kca_states(dimensions);
    try {
        std::vector<std::vector<double>> transition_matrix;
        json_obj.at("transition_matrix").get_to(transition_matrix);
        kca_states.setTransitionMatrix(transition_matrix);

        std::vector<std::vector<double>> transition_covariance;
        json_obj.at("transition_covariance").get_to(transition_covariance);
        kca_states.setTransitionCovariance(transition_covariance);

        std::vector<double> current_state_mean;
        json_obj.at("current_state_mean").get_to(current_state_mean);
        kca_states.setCurrentStateMean(current_state_mean);

        std::vector<std::vector<double>> current_state_covariance;
        json_obj.at("current_state_covariance")
            .get_to(current_state_covariance);
        kca_states.setCurrentStateCovariance(current_state_covariance);

        std::vector<std::vector<double>> observation_matrix;
        json_obj.at("observation_matrix").get_to(observation_matrix);
        kca_states.setObservationMatrix(observation_matrix);

        const double observation_offset =
            json_obj.at("observation_offset").template get<double>();
        kca_states.setObservationOffset(observation_offset);

        kca_states.setInitialized();
        return kca_states;
    } catch (const nlohmann::json::exception &exc) {
        throw json_parse_error(exc.what());
    }
}
