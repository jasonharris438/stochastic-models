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

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace {

// Upper bound for any deserialized filter dimension. Large enough for any
// real filter system while preventing hostile JSON from driving pathological
// matrix allocations.
constexpr std::int64_t kMaxFilterDimension{1024};

int getValidatedDimension(const nlohmann::json& json_obj, const char* key) {
  const nlohmann::json& field = json_obj.at(key);
  if (!field.is_number_integer()) {
    throw json_parse_error(
        "Dimension field '" + std::string{key} + "' must be a JSON integer."
    );
  }
  const std::int64_t value = field.template get<std::int64_t>();
  if (value < 1 || value > kMaxFilterDimension) {
    throw json_parse_error(
        "Dimension field '" + std::string{key} + "' must be in [1, " +
        std::to_string(kMaxFilterDimension) + "]; got " +
        std::to_string(value) + "."
    );
  }
  return static_cast<int>(value);
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

std::vector<std::vector<double>> getValidatedMatrix(
    const nlohmann::json& json_obj,
    const char* key,
    const std::size_t rows,
    const std::size_t columns
) {
  std::vector<std::vector<double>> matrix_as_vectors;
  json_obj.at(key).get_to(matrix_as_vectors);
  if (matrix_as_vectors.size() != rows) {
    throw json_parse_error(
        "Field '" + std::string{key} + "' must have " + std::to_string(rows) +
        " rows; got " + std::to_string(matrix_as_vectors.size()) + "."
    );
  }
  for (const std::vector<double>& matrix_row : matrix_as_vectors) {
    if (matrix_row.size() != columns) {
      throw json_parse_error(
          "Field '" + std::string{key} + "' must have " +
          std::to_string(columns) + " columns in every row; got a row of "
          "length " + std::to_string(matrix_row.size()) + "."
      );
    }
  }
  return matrix_as_vectors;
}

std::vector<double> getValidatedVector(
    const nlohmann::json& json_obj, const char* key, const std::size_t length
) {
  std::vector<double> vector_as_vector;
  json_obj.at(key).get_to(vector_as_vector);
  if (vector_as_vector.size() != length) {
    throw json_parse_error(
        "Field '" + std::string{key} + "' must have length " +
        std::to_string(length) + "; got " +
        std::to_string(vector_as_vector.size()) + "."
    );
  }
  return vector_as_vector;
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
  for (int i{0}; i < boost_matrix.size1(); i++) {
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

    FilterSystemDimensions dimensions;
    dimensions.state_mean_dimension =
        getValidatedDimension(json_obj, "state_mean_dimension");
    dimensions.state_covariance_rows =
        getValidatedDimension(json_obj, "state_covariance_rows");
    dimensions.state_covariance_columns =
        getValidatedDimension(json_obj, "state_covariance_columns");
    dimensions.observation_matrix_rows =
        getValidatedDimension(json_obj, "observation_matrix_rows");
    dimensions.observation_matrix_columns =
        getValidatedDimension(json_obj, "observation_matrix_columns");
    dimensions.observation_covariance_rows =
        getValidatedDimension(json_obj, "observation_covariance_rows");
    dimensions.observation_covariance_columns =
        getValidatedDimension(json_obj, "observation_covariance_columns");
    dimensions.observation_offset =
        getValidatedNumber(json_obj, "observation_offset");

    return dimensions;
  } catch (const nlohmann::json::exception& exc) {
    throw json_parse_error(exc.what());
  }
}
const std::string FilterSystemDimensionsJsonAdapter::serialize(
    const FilterSystemDimensions& dimensions
) const {
  nlohmann::json json_obj = {
      {"state_mean_dimension", dimensions.state_mean_dimension},
      {"state_covariance_rows", dimensions.state_covariance_rows},
      {"state_covariance_columns", dimensions.state_covariance_columns},
      {"observation_matrix_rows", dimensions.observation_matrix_rows},
      {"observation_matrix_columns", dimensions.observation_matrix_columns},
      {"observation_covariance_rows", dimensions.observation_covariance_rows},
      {"observation_covariance_columns",
       dimensions.observation_covariance_columns},
      {"observation_offset", dimensions.observation_offset}
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
  try {
    const nlohmann::json json_obj = nlohmann::json::parse(state);

    const std::size_t state_rows =
        static_cast<std::size_t>(dimensions.state_covariance_rows);
    const std::size_t state_columns =
        static_cast<std::size_t>(dimensions.state_covariance_columns);
    const std::size_t mean_length =
        static_cast<std::size_t>(dimensions.state_mean_dimension);
    const std::size_t observation_rows =
        static_cast<std::size_t>(dimensions.observation_matrix_rows);
    const std::size_t observation_columns =
        static_cast<std::size_t>(dimensions.observation_matrix_columns);

    KcaStates kca_states(dimensions);

    std::vector<std::vector<double>> transition_matrix = getValidatedMatrix(
        json_obj, "transition_matrix", state_rows, state_columns
    );
    kca_states.setTransitionMatrix(transition_matrix);

    std::vector<std::vector<double>> transition_covariance =
        getValidatedMatrix(
            json_obj, "transition_covariance", state_rows, state_columns
        );
    kca_states.setTransitionCovariance(transition_covariance);

    std::vector<double> current_state_mean =
        getValidatedVector(json_obj, "current_state_mean", mean_length);
    kca_states.setCurrentStateMean(current_state_mean);

    std::vector<std::vector<double>> current_state_covariance =
        getValidatedMatrix(
            json_obj, "current_state_covariance", state_rows, state_columns
        );
    kca_states.setCurrentStateCovariance(current_state_covariance);

    std::vector<std::vector<double>> observation_matrix = getValidatedMatrix(
        json_obj, "observation_matrix", observation_rows, observation_columns
    );
    kca_states.setObservationMatrix(observation_matrix);

    kca_states.setObservationOffset(
        getValidatedNumber(json_obj, "observation_offset")
    );

    kca_states.setInitialized();
    return kca_states;
  } catch (const nlohmann::json::exception& exc) {
    throw json_parse_error(exc.what());
  }
}
