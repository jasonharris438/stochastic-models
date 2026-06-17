// Include-hygiene guard for <stochastic_models/numeric_utils/helpers.h>.
//
// Context (issue #26): helpers.h used to include <gsl/gsl_integration.h>, which
// transitively pulls in <math.h>. In C++ that injects an `abs(double)` overload
// into the GLOBAL namespace. Because helpers.h is included by nearly every
// translation unit (it declares roundToDecimals), that leak silently changed
// how an unqualified `abs(double)` call resolved: with the leak it bound to the
// floating-point overload; without it, `abs` falls back to C's integer abs(int)
// from <cstdlib>, truncating the argument and returning 0 for any |x| < 1.
//
// We removed the GSL includes from helpers.h so it no longer leaks a global
// abs. This guard locks that hygiene in place. The detector below is evaluated
// against the include state produced by helpers.h ALONE (before <gtest/gtest.h>
// is pulled in, since gtest legitimately includes <cstdlib>/<cmath>). If
// helpers.h — or anything it includes — ever re-introduces an unqualified
// `abs(double)`, the static_assert fails and the test target FAILS TO COMPILE,
// turning a silent, hard-to-spot resolution change back into a loud build
// break.
//
// If this file stops compiling on the static_assert below, do NOT weaken the
// guard: find what re-added a transitive <math.h>/GSL include to helpers.h and
// move it into the .cpp that needs it (include-what-you-use).

#include "stochastic_models/numeric_utils/helpers.h"

// Unqualified call inside requires{} → ordinary + ADL lookup, evaluated as a
// soft boolean. True only if some `abs(double)` is visible in this TU, which
// (given only helpers.h is included so far) means helpers.h leaked one.
//
// The guard is deliberately `abs`-specific: that is the exact symbol whose
// integer/double overload resolution caused the issue-#26 bug. A broader probe
// (fabs, sqrt, ...) is intentionally avoided so the guard does not misfire if
// helpers.h ever legitimately needs <cmath> for some non-abs reason.
template <class T>
concept HelpersLeaksUnqualifiedAbs = requires(T x) { abs(x); };

static_assert(
    !HelpersLeaksUnqualifiedAbs<double>,
    "helpers.h leaks an unqualified abs(double) into its includers "
    "(include-hygiene regression; see issue #26). Move the offending "
    "transitive <math.h>/GSL include out of helpers.h into the .cpp that "
    "needs it."
);

#include <gtest/gtest.h>

// The real assertion this file makes is the COMPILE-TIME static_assert above;
// this case simply gives the guard a visible entry in the test report.
TEST(HelpersIncludeHygieneTest, HelpersDoesNotLeakGlobalAbs) {
  SUCCEED() << "helpers.h does not leak a global abs(double) (compile-time "
               "static_assert above enforces this).";
}
