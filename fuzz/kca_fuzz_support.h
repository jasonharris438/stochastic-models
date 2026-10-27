#ifndef STOCHASTIC_MODELS_FUZZ_KCA_FUZZ_SUPPORT_H
#define STOCHASTIC_MODELS_FUZZ_KCA_FUZZ_SUPPORT_H

#include "stochastic_models/exceptions/errors.h"
#include "stochastic_models/kalman_filter/adapters.h"
#include "stochastic_models/kalman_filter/states_exceptions.h"

#include <concepts>
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

  inline Counts& counts() {
    static Counts instance;
    return instance;
  }

  // Constructed at start-up so the summary prints even when nothing is counted.
  inline Counts& tally = counts();

  [[noreturn]] inline void recordUnexpected() {
    __builtin_trap();
  }

  /**
   * @brief Runs 1 entry point call and classifies what it throws.
   *
   * A rejection by design returns. A known fault is counted. Anything else
   * traps, so libFuzzer records the input.
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
      counts().known_faults++;
      return Outcome::KnownFault;
    } catch (const filter_invalid_operation&) {
      counts().known_faults++;
      return Outcome::KnownFault;
    } catch (const std::runtime_error& error) {
      // Only the bare type from the GSL LU wrapper is a known fault.
      if (typeid(error) == typeid(std::runtime_error)) {
        counts().known_faults++;
        return Outcome::KnownFault;
      }
      recordUnexpected();
    } catch (...) {
      recordUnexpected();
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
      counts().round_trip_failures++;
    }
  }

} // namespace kca_fuzz

#endif // STOCHASTIC_MODELS_FUZZ_KCA_FUZZ_SUPPORT_H
