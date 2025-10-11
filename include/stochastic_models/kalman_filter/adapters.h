#include "stochastic_models/kalman_filter/states.h"
#include "stochastic_models/kalman_filter/type_conversion.h"

#include <string>

/**
 * @brief Handles serialization and deserialization to and from JSON structure
 * for KcaStates objects.
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
  const KcaStates deserialize(
      const std::string& state, const FilterSystemDimensions& dimensions
  ) const;
};
/**
 * @brief Handles serialization and deserialization to and from JSON structure
 * for KcaParameters objects.
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
   * @brief Deserializes a JSON string to a FilterSystemDimensions object.
   *
   * Creates a readable JSON object from the `state` object and extracts
   * fields to initialize a FilterSystemDimensions object.
   */
  const FilterSystemDimensions deserialize(const std::string& state) const;
};
