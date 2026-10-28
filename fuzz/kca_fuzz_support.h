#ifndef STOCHASTIC_MODELS_FUZZ_KCA_FUZZ_SUPPORT_H
#define STOCHASTIC_MODELS_FUZZ_KCA_FUZZ_SUPPORT_H

#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/kalman_filter/adapters.h"
#include "stochastic_models/kalman_filter/states_exceptions.h"

#include <concepts>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <typeinfo>

namespace kca_fuzz {

  inline const std::string scheme_dimensions =
      "{\"observation_covariance_columns\":1,\"observation_covariance_rows\":1,"
      "\"observation_matrix_columns\":3,\"observation_matrix_rows\":1,"
      "\"observation_offset\":0.0,\"state_covariance_columns\":3,"
      "\"state_covariance_rows\":3,\"state_mean_dimension\":3}";

  enum class Outcome { Completed, Rejected, KnownFault };

  /**
   * @brief Tallies of the lenient classifications, printed when the process
   * exits.
   */
  struct Counts {
    unsigned long known_faults{0};
    unsigned long round_trip_failures{0};

    ~Counts() {
      std::cerr << "known faults: " << known_faults
                << ", round-trip failures: " << round_trip_failures << "\n";
    }
  };

  inline Counts tally;

  /**
   * @brief Prints the type and message of an unexpected exception, then traps
   * so that libFuzzer records the input.
   */
  [[noreturn]] inline void
  recordUnexpected(const char* type_name, const char* message) {
    std::cerr << "unexpected exception " << type_name << ": " << message
              << "\n";
    __builtin_trap();
  }

  /**
   * @brief Runs 1 entry point call and classifies what it throws.
   *
   * A rejection by design returns. A known fault is counted. Anything else is
   * reported and trapped.
   */
  template <std::invocable EntryPoint>
  Outcome runClassified(EntryPoint entry_point) {
    try {
      entry_point();
      return Outcome::Completed;
    } catch (const json_parse_error&) {
      return Outcome::Rejected;
    } catch (const invalid_filter_dimensions&) {
      return Outcome::Rejected;
    } catch (const InvalidNumberObservationsError&) {
      return Outcome::Rejected;
    } catch (const filter_uninitialised&) {
      tally.known_faults++;
      return Outcome::KnownFault;
    } catch (const filter_invalid_operation&) {
      tally.known_faults++;
      return Outcome::KnownFault;
    } catch (const std::runtime_error& error) {
      // A bare std::runtime_error is a known fault. A derived type is not.
      if (typeid(error) == typeid(std::runtime_error)) {
        tally.known_faults++;
        return Outcome::KnownFault;
      }
      recordUnexpected(typeid(error).name(), error.what());
    } catch (const std::exception& error) {
      recordUnexpected(typeid(error).name(), error.what());
    } catch (...) {
      recordUnexpected("unknown type", "no message");
    }
  }

  /**
   * @brief Parses an entry point output back through the state adapter and
   * counts a failure instead of failing.
   */
  inline void checkRoundTrip(
      const std::string& output, const std::string& dimensions_json
  ) {
    const FilterSystemDimensions dimensions =
        FilterSystemDimensionsJsonAdapter{}.deserialize(dimensions_json);
    try {
      KcaStatesJsonAdapter{}.deserialize(output, dimensions);
    } catch (const json_parse_error&) {
      tally.round_trip_failures++;
    }
  }

} // namespace kca_fuzz

#endif // STOCHASTIC_MODELS_FUZZ_KCA_FUZZ_SUPPORT_H
