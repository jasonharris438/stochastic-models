# Derivations

`expected_values.py` recomputes the expected values that the unit tests lock.
It uses the model formulas with numpy and scipy. It never calls the library.

It covers the Ornstein-Uhlenbeck, hitting time, trading level, value function,
Kalman filter and closed-form variance tests. It does not cover the General
Linear estimator outputs in `general_linear_likelihood_test.cpp`,
`general_linear_online_test.cpp` and `general_sde_model_test.cpp`. Those tests
lock the library's current output until the estimator is corrected.

The script writes `tests/support/expected_values.h`. The tests include that
header and name each value. Never edit the header by hand.

The script requires `uv`. Install it with:

    curl -LsSf https://astral.sh/uv/install.sh | sh

The dependencies are pinned in the script and in `expected_values.py.lock`.
Run from the repository root:

    uv run --locked --script tools/derivations/expected_values.py

That prints one line per value. To regenerate the header:

    uv run --locked --script tools/derivations/expected_values.py --header tests/support/expected_values.h

CI regenerates the header and fails when it differs from the committed file.
When a test value must change, change the formula or the inputs in the
script, regenerate the header, and commit both.
