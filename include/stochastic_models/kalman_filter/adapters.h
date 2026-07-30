#ifndef STOCHASTIC_MODELS_KALMAN_FILTER_ADAPTERS_H
#define STOCHASTIC_MODELS_KALMAN_FILTER_ADAPTERS_H

#include "stochastic_models/kalman_filter/states.h"
#include "stochastic_models/kalman_filter/type_conversion.h"

#include <string>

/**
 * @file
 * @brief JSON adapters used to serialize/deserialize Kalman filter state and
 * dimension objects.
 */

/**
 * @brief Handles serialization and deserialization to and from JSON for
 * `KcaStates` objects.
 */
class KcaStatesJsonAdapter {
private:
  const std::vector<std::vector<double>> copyBoostMatrixToVector(
      const boost::numeric::ublas::matrix<double>& boost_matrix
  ) const;
  const std::vector<double> copyBoostVectorToVector(
      const boost::numeric::ublas::vector<double>& boost_vector
  ) const;

public:
  /**
   * @brief Serializes a KcaStates object to a JSON string.
   *
   * @param kca_states The KcaStates object to serialize.
   * @return std::string The JSON string representation of the KcaStates
   * object.
   */
  const std::string serialize(const KcaStates& kca_states) const;
  /**
   * @brief Deserialize a JSON string into a `KcaStates` instance.
   *
   * @param state JSON string containing the state fields.
   * @param dimensions Dimensions object used to size the internal matrices.
   * @return KcaStates Reconstructed KCA state.
   * @throws json_parse_error If the JSON is unparseable, a field is missing
   *         or has the wrong type, or any array's shape does not match
   *         `dimensions`.
   */
  const KcaStates deserialize(
      const std::string& state, const FilterSystemDimensions& dimensions
  ) const;
};
/**
 * @brief Adapter to (de)serialize `FilterSystemDimensions` to/from JSON.
 */
class FilterSystemDimensionsJsonAdapter {
public:
  /**
   * @brief Serializes a FilterSystemDimensions object to a JSON string.
   *
   * @param dimensions The FilterSystemDimensions object to serialize.
   * @return std::string The JSON string representation of the
   * FilterSystemDimensions object.
   */
  const std::string serialize(const FilterSystemDimensions& dimensions) const;
  /**
   * @brief Parse a JSON string and return a `FilterSystemDimensions` object.
   *
   * The function expects the JSON to contain numeric fields matching the
   * struct layout.
   *
   * @param state JSON string to parse.
   * @return FilterSystemDimensions Parsed object.
   * @throws json_parse_error If the JSON is unparseable, a field is missing,
   *         a dimension is not an integer in [1, 1024], or the offset is not
   *         a number.
   */
  const FilterSystemDimensions deserialize(const std::string& state) const;
};

#endif // STOCHASTIC_MODELS_KALMAN_FILTER_ADAPTERS_H
