#ifndef STOCHASTIC_MODELS_TESTS_SUPPORT_CONVERSIONS_H
#define STOCHASTIC_MODELS_TESTS_SUPPORT_CONVERSIONS_H

#include <boost/numeric/ublas/vector.hpp>
#include <vector>

namespace test_support {

  inline std::vector<double>
  toStdVector(const boost::numeric::ublas::vector<double>& boost_vector) {
    return std::vector<double>(boost_vector.begin(), boost_vector.end());
  }

} // namespace test_support

#endif // STOCHASTIC_MODELS_TESTS_SUPPORT_CONVERSIONS_H
