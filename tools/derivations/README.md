# Derivations

`expected_values.py` computes the expected values in
`tests/support/expected_values.h` from the model formulas with numpy and scipy.
It never calls the library.

It does not cover the General Linear estimator outputs in
`general_linear_likelihood_test.cpp`, `general_linear_online_test.cpp` and
`general_sde_model_test.cpp`. Those tests lock the library's current output
until the estimator is corrected.

The tests include the header and name each value. Never edit the header by
hand.

The commands below need `uv` 0.12.7 or later. Install it from your package
manager, or download the release archive for your platform from
https://github.com/astral-sh/uv/releases and check it against its `.sha256`
file before you extract it. CI installs it that way.

The script needs Python 3.12 or later. CI runs it on Python 3.14. Without
`uv`, a virtual environment with `numpy==2.5.3` and `scipy==1.18.1` runs the
script too.

The dependencies are pinned in the script and in `expected_values.py.lock`.
Run from the repository root:

    uv run --locked --script tools/derivations/expected_values.py

That prints 1 line per value. To regenerate the header:

    uv run --locked --script tools/derivations/expected_values.py --header tests/support/expected_values.h

To check the committed header:

    uv run --locked --script tools/derivations/expected_values.py --check tests/support/expected_values.h

CI runs the check. It fails when a value differs from the script by more than
1e-10 relative, or when any other text differs. A difference in the last
digits of a value comes from the CPU that ran the script. Regenerate and
commit the header only when a formula or an input changes.
