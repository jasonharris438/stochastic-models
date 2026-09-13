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
#include <cstddef>
#include <cstdint>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

namespace {

  std::size_t
  getValidatedDimension(const nlohmann::json& json_obj, const char* key) {
    const nlohmann::json& field = json_obj.at(key);
    if (!field.is_number_integer()) {
      throw json_parse_error(
          "Dimension field '" + std::string{key} + "' must be a JSON integer."
      );
    }
    const std::int64_t value = field.template get<std::int64_t>();
    if (value < 0) {
      throw json_parse_error(
          "Dimension field '" + std::string{key} +
          "' must not be negative; got " + std::to_string(value) + "."
      );
    }
    return static_cast<std::size_t>(value);
  }

  double getValidatedNumber(const nlohmann::json& json_obj, const char* key) {
    const nlohmann::json& field = json_obj.at(key);
    if (!field.is_number()) {
      throw json_parse_error(
          "Field '" + std::string{key} + "' must be a JSON number."
      );
    }
    return field.template get<double>();
  }

} // namespace

const std::vector<std::vector<double>>
KcaStatesJsonAdapter::copyBoostMatrixToVector(
    const boost::numeric::ublas::matrix<double>& boost_matrix
) const {
  std::vector<std::vector<double>> result(
      boost_matrix.size1(), std::vector<double>(boost_matrix.size2())
  );

  // Copy the values from the matrix to the std::vector by row.
  for (std::size_t i{0}; i < boost_matrix.size1(); i++) {
    std::copy(
        row(boost_matrix, i).begin(), row(boost_matrix, i).end(),
        result.at(i).begin()
    );
  }
  return result;
}
const std::vector<double> KcaStatesJsonAdapter::copyBoostVectorToVector(
    const boost::numeric::ublas::vector<double>& boost_vector
) const {
  std::vector<double> std_vector(boost_vector.size());
  std::copy(boost_vector.begin(), boost_vector.end(), std_vector.begin());
  return std_vector;
}
const FilterSystemDimensions
FilterSystemDimensionsJsonAdapter::deserialize(const std::string& state) const {
  try {
    const nlohmann::json json_obj = nlohmann::json::parse(state);

    const std::size_t state_mean_dimension =
        getValidatedDimension(json_obj, "state_mean_dimension");
    const std::size_t state_covariance_rows =
        getValidatedDimension(json_obj, "state_covariance_rows");
    const std::size_t state_covariance_columns =
        getValidatedDimension(json_obj, "state_covariance_columns");
    const std::size_t observation_matrix_rows =
        getValidatedDimension(json_obj, "observation_matrix_rows");
    const std::size_t observation_matrix_columns =
        getValidatedDimension(json_obj, "observation_matrix_columns");
    const std::size_t observation_covariance_rows =
        getValidatedDimension(json_obj, "observation_covariance_rows");
    const std::size_t observation_covariance_columns =
        getValidatedDimension(json_obj, "observation_covariance_columns");
    const double observation_offset =
        getValidatedNumber(json_obj, "observation_offset");

    return FilterSystemDimensions(
        state_mean_dimension, state_covariance_rows, state_covariance_columns,
        observation_matrix_rows, observation_matrix_columns,
        observation_covariance_rows, observation_covariance_columns,
        observation_offset
    );
  } catch (const nlohmann::json::exception& exc) {
    throw json_parse_error(exc.what());
  } catch (const invalid_filter_dimensions& exc) {
    throw json_parse_error(exc.what());
  }
}
const std::string FilterSystemDimensionsJsonAdapter::serialize(
    const FilterSystemDimensions& dimensions
) const {
  nlohmann::json json_obj = {
      {"state_mean_dimension", dimensions.getStateMeanDimension()},
      {"state_covariance_rows", dimensions.getStateCovarianceRows()},
      {"state_covariance_columns", dimensions.getStateCovarianceColumns()},
      {"observation_matrix_rows", dimensions.getObservationMatrixRows()},
      {"observation_matrix_columns", dimensions.getObservationMatrixColumns()},
      {"observation_covariance_rows",
       dimensions.getObservationCovarianceRows()},
      {"observation_covariance_columns",
       dimensions.getObservationCovarianceColumns()},
      {"observation_offset", dimensions.getObservationOffset()}
  };
  return json_obj.dump();
}
const std::string
KcaStatesJsonAdapter::serialize(const KcaStates& kca_states) const {
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
    const std::string& state, const FilterSystemDimensions& dimensions
) const {
  KcaStates kca_states(dimensions);
  try {
    const nlohmann::json json_obj = nlohmann::json::parse(state);

    kca_states.setTransitionMatrix(
        json_obj.at("transition_matrix").get<std::vector<std::vector<double>>>()
    );
    kca_states.setTransitionCovariance(
        json_obj.at("transition_covariance")
            .get<std::vector<std::vector<double>>>()
    );
    kca_states.setCurrentStateMean(
        json_obj.at("current_state_mean").get<std::vector<double>>()
    );
    kca_states.setCurrentStateCovariance(
        json_obj.at("current_state_covariance")
            .get<std::vector<std::vector<double>>>()
    );
    kca_states.setObservationMatrix(
        json_obj.at("observation_matrix")
            .get<std::vector<std::vector<double>>>()
    );

    kca_states.setObservationOffset(
        getValidatedNumber(json_obj, "observation_offset")
    );

    kca_states.setInitialized();
    return kca_states;
  } catch (const nlohmann::json::exception& exc) {
    throw json_parse_error(exc.what());
  } catch (const filter_shape_mismatch& exc) {
    throw json_parse_error(exc.what());
  }
}
