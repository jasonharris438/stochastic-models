#include "stochastic_models/numeric_utils/linalg.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <gsl/gsl_errno.h>
#include <gtest/gtest.h>
#include <stdexcept>

namespace {
  // Reads the currently-installed GSL handler without lasting side effects:
  // gsl_set_error_handler_off() installs a no-op and returns the previous
  // handler, which we then re-install. The GSL default handler reads back as
  // nullptr.
  gsl_error_handler_t* current_gsl_handler() {
    gsl_error_handler_t* h = gsl_set_error_handler_off();
    gsl_set_error_handler(h);
    return h;
  }
} // namespace

/**
 * @test When a guarded GSL routine throws, the previously-installed error
 * handler must be restored (RAII guard). A singular matrix forces
 * BoostMatrixInverter::invertMatrix to throw mid-routine.
 */
TEST(GslHandlerGuardTest, HandlerRestoredAfterThrow) {
  namespace ublas = boost::numeric::ublas;

  gsl_error_handler_t* before = current_gsl_handler();

  ublas::matrix<double> singular(2, 2);
  singular(0, 0) = 1.0;
  singular(0, 1) = 2.0;
  singular(1, 0) = 2.0;
  singular(1, 1) = 4.0; // rows linearly dependent

  BoostMatrixInverter inverter;
  EXPECT_THROW(inverter.invertMatrix(singular), std::runtime_error);

  gsl_error_handler_t* after = current_gsl_handler();
  EXPECT_EQ(
      after, before
  ) << "GSL error handler was not restored after a throw (handler leaked).";
}
