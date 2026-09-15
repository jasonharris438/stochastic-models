# Derivations

`expected_values.py` recomputes every numeric value that the unit tests lock.
It uses the model formulas with numpy and scipy. It never calls the library.

The script requires `uv`. It installs numpy and scipy in an isolated
environment for each run. Install `uv` with:

    curl -LsSf https://astral.sh/uv/install.sh | sh

Run the script from the repository root:

    uv run --with numpy --with scipy tools/derivations/expected_values.py

Each line prints a test name and the value at full double precision. When a
test value changes, run the script and compare before you edit the test.
