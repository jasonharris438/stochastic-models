# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy==2.5.3", "scipy==1.18.1"]
# ///
"""Independent replica of the expected values in tests/support/expected_values.h.

Every value is computed from the model formulas with scipy, not from the
library. Quadrature uses QAWS for the algebraic endpoint singularity and
brentq for the roots. Run with:

    uv run --script tools/derivations/expected_values.py

Without arguments the script prints 1 `name value` line per expected value,
keyed by the gtest identifier that locks it. With `--header PATH` it writes the
same values as a C++ header of compile time constants. With `--check PATH` it
compares the header at PATH with the same values and exits 1 when a value is
not within 1e-10 relative or when any other text differs.
"""

import argparse
import math
import re
import sys
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm

ExpectedValue = float | int | list[float] | list[list[float]]


@dataclass(frozen=True)
class OrnsteinUhlenbeck:
    """Parameters of dX = alpha (mu - X) dt + sigma dW."""

    mu: float
    alpha: float
    sigma: float

    @property
    def stationary_deviation(self) -> float:
        """Standard deviation of the stationary law, sigma / sqrt(2 alpha)."""
        return self.sigma / np.sqrt(2.0 * self.alpha)

    @property
    def scaled_drift(self) -> float:
        """The factor sqrt(2 alpha / sigma^2) in the trading kernels."""
        return np.sqrt(2.0 * self.alpha) / self.sigma

    def conditional_variance(self, step: float) -> float:
        """Variance over 1 step, sigma^2 (1 - exp(-2 alpha t)) / (2 alpha)."""
        if abs(self.alpha) < 1e-12:
            return self.sigma**2 * step
        return float(self.sigma**2 * -np.expm1(-2.0 * self.alpha * step) / (2.0 * self.alpha))


@dataclass(frozen=True)
class GeneralLinear:
    """Parameters of dX = mu X dt + sigma dW."""

    mu: float
    sigma: float

    @property
    def unconditional_variance(self) -> float:
        """Stationary variance sigma^2 / (-2 mu), defined only for mu < 0."""
        return self.sigma**2 / (-2.0 * self.mu)

    def conditional_variance(self, step: float) -> float:
        """Variance over 1 step, sigma^2 (exp(2 mu t) - 1) / (2 mu)."""
        if abs(self.mu) < 1e-12:
            return self.sigma**2 * step
        return float(self.sigma**2 * np.expm1(2.0 * self.mu * step) / (2.0 * self.mu))


@dataclass(frozen=True)
class TradingCosts:
    """Discount rate and transaction cost of the trading problem."""

    rate: float
    cost: float


@dataclass(frozen=True)
class OrnsteinUhlenbeckSums:
    """Sufficient statistics over the lag and lead pairs of a series."""

    lead_sum: float
    lag_sum: float
    lead_sum_squared: float
    lag_sum_squared: float
    lead_lag_sum_product: float
    n_obs: int


def hitting_time_core(point: float, *, model: OrnsteinUhlenbeck) -> float:
    """Unnormalised hitting time kernel exp(x alpha (x - 2 mu) / sigma^2)."""
    return float(np.exp(x_exponent(point, model=model)))


def hitting_time_density(
    point: float, *, model: OrnsteinUhlenbeck, first: float, second: float
) -> float:
    """Ratio of the kernel integral over [second, point] to that over [second, first].

    The exponent is shifted by its value at `second` so that both integrals
    stay in floating point range.
    """
    shift = x_exponent(second, model=model)

    def shifted_kernel(point: float) -> float:
        return float(np.exp(x_exponent(point, model=model) - shift))

    numerator, _ = quad(shifted_kernel, second, point, epsabs=0.0, epsrel=1e-12)
    denominator, _ = quad(shifted_kernel, second, first, epsabs=0.0, epsrel=1e-12)
    return numerator / denominator


def x_exponent(point: float, *, model: OrnsteinUhlenbeck) -> float:
    """Exponent of the hitting time kernel at 1 point."""
    return point * model.alpha * (point - 2.0 * model.mu) / model.sigma**2


def kernel_integral(
    point: float,
    *,
    model: OrnsteinUhlenbeck,
    rate: float,
    mean_sign: float,
    power_shift: int = 0,
    lower: float = 0.0,
) -> float:
    """Integral of u^(r/alpha - 1 + shift) exp(sign k (x - mu) u - u^2 / 2) over [lower, inf).

    With sign = +1 this is the Leung-Li F function, with sign = -1 the G
    function. A power shift of 1 gives the derivative integrand. The
    algebraic singularity at zero is handled by QAWS on [0, 1].
    """
    power = rate / model.alpha - 1.0 + power_shift
    drift = mean_sign * model.scaled_drift * (point - model.mu)

    def smooth_part(integration_variable: float) -> float:
        return float(
            np.exp(drift * integration_variable - 0.5 * integration_variable * integration_variable)
        )

    def full_integrand(integration_variable: float) -> float:
        return integration_variable**power * smooth_part(integration_variable)

    if lower > 0.0:
        tail, _ = quad(full_integrand, lower, np.inf, epsabs=0.0, epsrel=1e-12, limit=200)
        return tail
    head, _ = quad(smooth_part, 0.0, 1.0, weight="alg", wvar=(power, 0.0), epsabs=0.0, epsrel=1e-12)
    tail, _ = quad(full_integrand, 1.0, np.inf, epsabs=0.0, epsrel=1e-12, limit=200)
    return head + tail


def function_f(point: float, *, model: OrnsteinUhlenbeck, rate: float) -> float:
    """Leung-Li F(x)."""
    return kernel_integral(point, model=model, rate=rate, mean_sign=1.0)


def function_g(point: float, *, model: OrnsteinUhlenbeck, rate: float) -> float:
    """Leung-Li G(x)."""
    return kernel_integral(point, model=model, rate=rate, mean_sign=-1.0)


def derivative_f(point: float, *, model: OrnsteinUhlenbeck, rate: float) -> float:
    """F'(x) by differentiation under the integral sign."""
    return model.scaled_drift * kernel_integral(
        point, model=model, rate=rate, mean_sign=1.0, power_shift=1
    )


def derivative_g(point: float, *, model: OrnsteinUhlenbeck, rate: float) -> float:
    """G'(x) by differentiation under the integral sign."""
    return -model.scaled_drift * kernel_integral(
        point, model=model, rate=rate, mean_sign=-1.0, power_shift=1
    )


def level_l_star(*, model: OrnsteinUhlenbeck, costs: TradingCosts) -> float:
    """The level L* = (alpha mu + r c) / (r + alpha)."""
    return (model.alpha * model.mu + costs.rate * costs.cost) / (costs.rate + model.alpha)


def exit_residual(level: float, *, model: OrnsteinUhlenbeck, costs: TradingCosts) -> float:
    """Smooth pasting residual F(b) - (b - c) F'(b) of the exit problem."""
    return function_f(level, model=model, rate=costs.rate) - (level - costs.cost) * derivative_f(
        level, model=model, rate=costs.rate
    )


def exit_residual_exponential(
    level: float, *, model: OrnsteinUhlenbeck, costs: TradingCosts
) -> float:
    """Smooth pasting residual e^b F(b) - (e^b - c) F'(b) of the exponential exit problem."""
    return np.exp(level) * function_f(level, model=model, rate=costs.rate) - (
        np.exp(level) - costs.cost
    ) * derivative_f(level, model=model, rate=costs.rate)


def exit_residual_stop_loss(
    level: float, *, model: OrnsteinUhlenbeck, costs: TradingCosts, stop_loss: float
) -> float:
    """Smooth pasting residual of the exit problem with a stop loss.

    The value function C F + D G matches L - c at the stop loss and b - c at
    the exit. The residual is the smooth pasting condition at b, multiplied
    by the common denominator, in the sign convention the library uses.
    """
    rate = costs.rate
    level_less_cost = level - costs.cost
    stop_less_cost = stop_loss - costs.cost
    g_stop = function_g(stop_loss, model=model, rate=rate)
    g_level = function_g(level, model=model, rate=rate)
    f_stop = function_f(stop_loss, model=model, rate=rate)
    f_level = function_f(level, model=model, rate=rate)
    return (
        (stop_less_cost * g_level - level_less_cost * g_stop)
        * derivative_f(level, model=model, rate=rate)
        + (level_less_cost * f_stop - stop_less_cost * f_level)
        * derivative_g(level, model=model, rate=rate)
        - (g_level * f_stop - g_stop * f_level)
    )


def value_function(
    point: float, *, model: OrnsteinUhlenbeck, costs: TradingCosts, exit_level: float
) -> float:
    """Value function (b - c) F(x) / F(b) below the exit, x - c above it."""
    if point >= exit_level:
        return point - costs.cost
    rate = costs.rate
    return (
        (exit_level - costs.cost)
        * function_f(point, model=model, rate=rate)
        / function_f(exit_level, model=model, rate=rate)
    )


def derivative_value_function(
    point: float, *, model: OrnsteinUhlenbeck, costs: TradingCosts, exit_level: float
) -> float:
    """V'(x) = (b - c) F'(x) / F(b) below the exit, 1 above it."""
    if point >= exit_level:
        return 1.0
    rate = costs.rate
    return (
        (exit_level - costs.cost)
        * derivative_f(point, model=model, rate=rate)
        / function_f(exit_level, model=model, rate=rate)
    )


def stop_loss_weights(
    *, model: OrnsteinUhlenbeck, costs: TradingCosts, exit_level: float, stop_loss: float
) -> tuple[float, float]:
    """Weights C and D of the stop loss value function C F + D G."""
    rate = costs.rate
    level_less_cost = exit_level - costs.cost
    stop_less_cost = stop_loss - costs.cost
    f_level = function_f(exit_level, model=model, rate=rate)
    f_stop = function_f(stop_loss, model=model, rate=rate)
    g_level = function_g(exit_level, model=model, rate=rate)
    g_stop = function_g(stop_loss, model=model, rate=rate)
    denominator = f_level * g_stop - f_stop * g_level
    weight_f = (level_less_cost * g_stop - stop_less_cost * g_level) / denominator
    weight_g = (stop_less_cost * f_level - level_less_cost * f_stop) / denominator
    return weight_f, weight_g


def value_function_stop_loss(
    point: float,
    *,
    model: OrnsteinUhlenbeck,
    costs: TradingCosts,
    exit_level: float,
    stop_loss: float,
) -> float:
    """Value function C F(x) + D G(x) between the stop loss and the exit, x - c outside."""
    if not (exit_level > point > stop_loss):
        return point - costs.cost
    weight_f, weight_g = stop_loss_weights(
        model=model, costs=costs, exit_level=exit_level, stop_loss=stop_loss
    )
    rate = costs.rate
    return weight_f * function_f(point, model=model, rate=rate) + weight_g * function_g(
        point, model=model, rate=rate
    )


def derivative_value_function_stop_loss(
    point: float,
    *,
    model: OrnsteinUhlenbeck,
    costs: TradingCosts,
    exit_level: float,
    stop_loss: float,
) -> float:
    """V'(x) = C F'(x) + D G'(x) from the stop loss up to the exit, 1 outside.

    At the stop loss itself the slope is the limit from inside the
    continuation region, where the payoff kink makes the left slope differ.
    """
    if not (exit_level > point >= stop_loss):
        return 1.0
    weight_f, weight_g = stop_loss_weights(
        model=model, costs=costs, exit_level=exit_level, stop_loss=stop_loss
    )
    rate = costs.rate
    return weight_f * derivative_f(point, model=model, rate=rate) + weight_g * derivative_g(
        point, model=model, rate=rate
    )


def value_function_exponential(
    point: float, *, model: OrnsteinUhlenbeck, costs: TradingCosts, exit_level: float
) -> float:
    """Value function (e^b - c) F(x) / F(b) below the exit, e^x - c above it."""
    if point >= exit_level:
        return float(np.exp(point)) - costs.cost
    rate = costs.rate
    return (
        (float(np.exp(exit_level)) - costs.cost)
        * function_f(point, model=model, rate=rate)
        / function_f(exit_level, model=model, rate=rate)
    )


def derivative_value_function_exponential(
    point: float, *, model: OrnsteinUhlenbeck, costs: TradingCosts, exit_level: float
) -> float:
    """V'(x) = (e^b - c) F'(x) / F(b) below the exit, e^x above it."""
    if point >= exit_level:
        return float(np.exp(point))
    rate = costs.rate
    return (
        (float(np.exp(exit_level)) - costs.cost)
        * derivative_f(point, model=model, rate=rate)
        / function_f(exit_level, model=model, rate=rate)
    )


def entry_residual(
    level: float, *, model: OrnsteinUhlenbeck, costs: TradingCosts, exit_level: float
) -> float:
    """Entry residual G(d) (V'(d) - 1) - G'(d) (V(d) - d - c) with no stop loss."""
    rate = costs.rate
    value = value_function(level, model=model, costs=costs, exit_level=exit_level)
    slope = derivative_value_function(level, model=model, costs=costs, exit_level=exit_level)
    return function_g(level, model=model, rate=rate) * (slope - 1.0) - derivative_g(
        level, model=model, rate=rate
    ) * (value - level - costs.cost)


def entry_residual_stop_loss(
    level: float,
    *,
    model: OrnsteinUhlenbeck,
    costs: TradingCosts,
    exit_level: float,
    stop_loss: float,
) -> float:
    """Entry residual G(d) (V'(d) - 1) - G'(d) (V(d) - d - c) with a stop loss."""
    rate = costs.rate
    value = value_function_stop_loss(
        level, model=model, costs=costs, exit_level=exit_level, stop_loss=stop_loss
    )
    slope = derivative_value_function_stop_loss(
        level, model=model, costs=costs, exit_level=exit_level, stop_loss=stop_loss
    )
    return function_g(level, model=model, rate=rate) * (slope - 1.0) - derivative_g(
        level, model=model, rate=rate
    ) * (value - level - costs.cost)


def lower_entry_residual_stop_loss(
    level: float,
    *,
    model: OrnsteinUhlenbeck,
    costs: TradingCosts,
    exit_level: float,
    stop_loss: float,
) -> float:
    """Lower entry residual F(a) (V'(a) - 1) - F'(a) (V(a) - a - c) with a stop loss."""
    rate = costs.rate
    value = value_function_stop_loss(
        level, model=model, costs=costs, exit_level=exit_level, stop_loss=stop_loss
    )
    slope = derivative_value_function_stop_loss(
        level, model=model, costs=costs, exit_level=exit_level, stop_loss=stop_loss
    )
    return function_f(level, model=model, rate=rate) * (slope - 1.0) - derivative_f(
        level, model=model, rate=rate
    ) * (value - level - costs.cost)


def entry_residual_exponential(
    level: float, *, model: OrnsteinUhlenbeck, costs: TradingCosts, exit_level: float
) -> float:
    """Entry residual G(d) (V'(d) - e^d) - G'(d) (V(d) - e^d - c) of the exponential model."""
    rate = costs.rate
    payoff = float(np.exp(level))
    value = value_function_exponential(level, model=model, costs=costs, exit_level=exit_level)
    slope = derivative_value_function_exponential(
        level, model=model, costs=costs, exit_level=exit_level
    )
    return function_g(level, model=model, rate=rate) * (slope - payoff) - derivative_g(
        level, model=model, rate=rate
    ) * (value - payoff - costs.cost)


def solve_level(residual: Callable[[float], float], *, lower: float, upper: float) -> float:
    """Root of the residual inside the library's solver bracket."""
    return float(brentq(residual, lower, upper, xtol=1e-14, rtol=1e-14, maxiter=500))


def solve_exit_level(
    residual: Callable[[float], float], *, model: OrnsteinUhlenbeck, costs: TradingCosts
) -> float:
    """Root of the exit residual inside the library's exit bracket."""
    return solve_level(
        residual,
        lower=max(level_l_star(model=model, costs=costs), costs.cost),
        upper=model.mu + 4.0 * model.stationary_deviation,
    )


def ornstein_uhlenbeck_sums(series: list[float]) -> OrnsteinUhlenbeckSums:
    """Sufficient statistics of a series with unit steps."""
    lag = np.asarray(series[:-1])
    lead = np.asarray(series[1:])
    return OrnsteinUhlenbeckSums(
        lead_sum=float(lead.sum()),
        lag_sum=float(lag.sum()),
        lead_sum_squared=float((lead**2).sum()),
        lag_sum_squared=float((lag**2).sum()),
        lead_lag_sum_product=float((lead * lag).sum()),
        n_obs=len(series),
    )


def update_sums(
    sums: OrnsteinUhlenbeckSums, *, new_observation: float, last_observation: float
) -> OrnsteinUhlenbeckSums:
    """Sufficient statistics after 1 more lag and lead pair."""
    return OrnsteinUhlenbeckSums(
        lead_sum=sums.lead_sum + new_observation,
        lag_sum=sums.lag_sum + last_observation,
        lead_sum_squared=sums.lead_sum_squared + new_observation**2,
        lag_sum_squared=sums.lag_sum_squared + last_observation**2,
        lead_lag_sum_product=sums.lead_lag_sum_product + new_observation * last_observation,
        n_obs=sums.n_obs + 1,
    )


def ornstein_uhlenbeck_estimate(sums: OrnsteinUhlenbeckSums) -> OrnsteinUhlenbeck:
    """Maximum likelihood estimate from the sufficient statistics, unit steps."""
    pairs = sums.n_obs - 1
    mu = (sums.lead_sum * sums.lag_sum_squared - sums.lag_sum * sums.lead_lag_sum_product) / (
        pairs * (sums.lag_sum_squared - sums.lead_lag_sum_product)
        - (sums.lag_sum**2 - sums.lag_sum * sums.lead_sum)
    )
    alpha = np.log(sums.lag_sum_squared - 2.0 * mu * sums.lag_sum + pairs * mu**2) - np.log(
        sums.lead_lag_sum_product - mu * (sums.lag_sum + sums.lead_sum) + pairs * mu**2
    )
    decay = np.exp(-alpha)
    transition_variance = (
        sums.lead_sum_squared
        - 2.0 * decay * sums.lead_lag_sum_product
        + decay**2 * sums.lag_sum_squared
        - 2.0 * mu * (1.0 - decay) * (sums.lead_sum - decay * sums.lag_sum)
        + pairs * mu**2 * (1.0 - decay) ** 2
    ) / pairs
    sigma = np.sqrt(transition_variance * 2.0 * alpha / (1.0 - decay**2))
    return OrnsteinUhlenbeck(mu=float(mu), alpha=float(alpha), sigma=float(sigma))


def vector_values(values: np.ndarray) -> list[float]:
    """Plain float list of a numpy vector."""
    return [float(value) for value in values]


def matrix_values(values: np.ndarray) -> list[list[float]]:
    """Plain float list of lists of a numpy matrix."""
    return [[float(value) for value in row] for row in values]


@dataclass(frozen=True)
class KcaState:
    """Kinetic components state: 3-state constant-acceleration Kalman filter."""

    transition_matrix: np.ndarray
    transition_covariance: np.ndarray
    current_state_mean: np.ndarray
    current_state_covariance: np.ndarray
    observation_matrix: np.ndarray
    observation_offset: float = 0.0


KCA_SERIES = [
    10.51255,
    10.51985,
    10.52405,
    10.4656,
    10.47,
    10.5403,
    10.4425,
    10.3087,
    10.1994,
    10.1839,
    10.24645,
    10.1795,
    10.21715,
    10.14995,
    10.194,
    10.22505,
    10.27325,
    10.25095,
    10.30575,
    10.27645,
]


@dataclass(frozen=True)
class KcaStep:
    """Intermediate quantities and result of 1 predict and correct step."""

    predicted_state_mean: np.ndarray
    predicted_state_covariance: np.ndarray
    predicted_observation_mean: np.ndarray
    predicted_observation_covariance: np.ndarray
    kalman_gain: np.ndarray
    state: KcaState


def kca_initial_state(series: list[float], *, step: float, noise: float) -> KcaState:
    """Initial state: last observation as position, zero covariance, white-noise jerk model."""
    transition = np.array([[1.0, step, step**2 / 2.0], [0.0, 1.0, step], [0.0, 0.0, 1.0]])
    covariance = noise * np.array(
        [
            [step**5 / 20.0, step**4 / 8.0, step**3 / 6.0],
            [step**4 / 8.0, step**3 / 3.0, step**2 / 2.0],
            [step**3 / 6.0, step**2 / 2.0, step],
        ]
    )
    return KcaState(
        transition_matrix=transition,
        transition_covariance=covariance,
        current_state_mean=np.array([series[-1], 0.0, 0.0]),
        current_state_covariance=np.zeros((3, 3)),
        observation_matrix=np.array([[1.0, 0.0, 0.0]]),
    )


def kca_prior_state() -> KcaState:
    """The prior state that the update tests set by hand."""
    return KcaState(
        transition_matrix=np.array(
            [[1.0011961162353782, 1.0, 0.5], [0.0, 1.0, 1.0], [0.0, 0.0, 1.0]]
        ),
        transition_covariance=np.diag([0.12695229227341848, 0.001, 0.001]),
        current_state_mean=np.array([10.288741828687053, 0.0, 0.0]),
        current_state_covariance=np.zeros((3, 3)),
        observation_matrix=np.array([[1.0, 0.0, 0.0]]),
    )


def kca_step(state: KcaState, *, observation: float, innovation_sigma: float) -> KcaStep:
    """1 predict and correct step of the Kalman filter."""
    transition = state.transition_matrix
    predicted_mean = transition @ state.current_state_mean
    predicted_covariance = (
        transition @ state.current_state_covariance @ transition.T + state.transition_covariance
    )
    observation_matrix = state.observation_matrix
    predicted_observation_mean = observation_matrix @ predicted_mean + state.observation_offset
    innovation_covariance = (
        observation_matrix @ predicted_covariance @ observation_matrix.T + innovation_sigma**2
    )
    gain = predicted_covariance @ observation_matrix.T @ np.linalg.inv(innovation_covariance)
    innovation = observation - predicted_observation_mean[0]
    corrected_mean = predicted_mean + gain[:, 0] * innovation
    corrected_covariance = predicted_covariance - gain @ (observation_matrix @ predicted_covariance)
    return KcaStep(
        predicted_state_mean=predicted_mean,
        predicted_state_covariance=predicted_covariance,
        predicted_observation_mean=predicted_observation_mean,
        predicted_observation_covariance=innovation_covariance,
        kalman_gain=gain,
        state=KcaState(
            transition_matrix=transition,
            transition_covariance=state.transition_covariance,
            current_state_mean=corrected_mean,
            current_state_covariance=corrected_covariance,
            observation_matrix=observation_matrix,
            observation_offset=state.observation_offset,
        ),
    )


def kca_update(state: KcaState, *, observation: float, innovation_sigma: float) -> KcaState:
    """State after 1 predict and correct step of the Kalman filter."""
    return kca_step(state, observation=observation, innovation_sigma=innovation_sigma).state


def kca_state_values(prefix: str, state: KcaState, fields: list[str]) -> dict[str, ExpectedValue]:
    """Expected value entries for the named fields of 1 KCA state."""
    values: dict[str, ExpectedValue] = {}
    for field in fields:
        member = getattr(state, field)
        if isinstance(member, float):
            values[f"{prefix}.{field}"] = member
        elif member.ndim == 1:
            values[f"{prefix}.{field}"] = vector_values(member)
        else:
            values[f"{prefix}.{field}"] = matrix_values(member)
    return values


def kca_expected_values() -> dict[str, ExpectedValue]:
    """Expected state fields of the KCA entry point, filter state and filter update tests."""
    initial = kca_initial_state(KCA_SERIES, step=1.0, noise=0.001)
    values = kca_state_values(
        "KcaTest.getInitializedKcaStateTest",
        initial,
        [
            "transition_matrix",
            "transition_covariance",
            "current_state_mean",
            "current_state_covariance",
            "observation_matrix",
            "observation_offset",
        ],
    )
    posterior = kca_update(kca_prior_state(), observation=10.3, innovation_sigma=0.1)
    values.update(
        kca_state_values(
            "KcaTest.getUpdatedKcaStateTest",
            posterior,
            [
                "current_state_mean",
                "current_state_covariance",
                "transition_matrix",
                "transition_covariance",
                "observation_matrix",
                "observation_offset",
            ],
        )
    )
    values["KalmanFilterStateTest.KcaStatessetInitialStateTest.transition_covariance"] = (
        matrix_values(initial.transition_covariance)
    )
    values["KalmanFilterUpdateTest.KineticComponentsupdatePosteriorsTest.current_state_mean"] = (
        vector_values(posterior.current_state_mean)
    )
    values["KalmanFilterUpdateTest.KineticComponentsFullFilterTest.current_state_mean"] = (
        vector_values(
            kca_update(initial, observation=10.3, innovation_sigma=0.1).current_state_mean
        )
    )
    return values


def filter_states_expected_values() -> dict[str, ExpectedValue]:
    """Expected vectors and matrices of the Kalman filter state unit tests."""
    values: dict[str, ExpectedValue] = {}
    observation_matrix = np.array([[1.0, 0.0, 0.0]])
    state_covariance = np.diag([0.013744, 0.001, 0.001])
    zero_covariance = np.zeros((3, 3))

    transition_matrix = np.array([[1.000295, 1.0, 0.5], [0.0, 1.0, 1.0], [0.0, 0.0, 1.0]])
    values["KalmanFilterTest.PredictedStateCalculateCovarianceTest.predicted_state_covariance"] = (
        matrix_values(transition_matrix @ zero_covariance @ transition_matrix.T + state_covariance)
    )
    values["KalmanFilterTest.PredictedStateCalculateMeanTest.predicted_state_mean"] = vector_values(
        transition_matrix @ np.array([1.330593, 0.0, 0.0])
    )

    predicted_state_mean = np.array([1.330986, 0.0, 0.0])
    values["KalmanFilterTest.PredictedObservationCalculateMeanTest.predicted_observation_mean"] = (
        vector_values(observation_matrix @ predicted_state_mean + 0.0)
    )
    values[
        "KalmanFilterTest.PredictedObservationCalculateCovarianceTest"
        ".predicted_observation_covariance"
    ] = matrix_values(observation_matrix @ state_covariance @ observation_matrix.T + 0.00687526**2)
    values["KalmanFilterTest.PredictedObservationCalculateKalmanGainTest.kalman_gain"] = (
        matrix_values(
            state_covariance @ observation_matrix.T @ np.linalg.inv(np.array([[0.0137917]]))
        )
    )

    gain_column = np.array([[0.99657263], [0.0], [0.0]])
    values["KalmanFilterTest.CurrentStateCalculateMeanTest.current_state_mean"] = vector_values(
        predicted_state_mean + gain_column[:, 0] * -0.02018567
    )
    values["KalmanFilterTest.CurrentStateCalculateCovarianceTest.current_state_covariance"] = (
        matrix_values(state_covariance - gain_column @ (observation_matrix @ state_covariance))
    )

    step = kca_step(kca_prior_state(), observation=10.3, innovation_sigma=0.1)
    values["KalmanFilterStateTest.KcaStatesupdatePredictedStateTest.predicted_state_mean"] = (
        vector_values(step.predicted_state_mean)
    )
    values["KalmanFilterStateTest.KcaStatesupdatePredictedStateTest.predicted_state_covariance"] = (
        matrix_values(step.predicted_state_covariance)
    )
    values["KalmanFilterStateTest.KcaStatesupdateCurrentStateTest.predicted_observation_mean"] = (
        vector_values(step.predicted_observation_mean)
    )
    values[
        "KalmanFilterStateTest.KcaStatesupdateCurrentStateTest.predicted_observation_covariance"
    ] = matrix_values(step.predicted_observation_covariance)
    values["KalmanFilterStateTest.KcaStatesupdateCurrentStateTest.current_state_mean"] = (
        vector_values(step.state.current_state_mean)
    )
    values["KalmanFilterStateTest.KcaStatesupdateCurrentStateTest.current_state_covariance"] = (
        matrix_values(step.state.current_state_covariance)
    )
    return values


def expected_values() -> dict[str, ExpectedValue]:
    """Every expected value, keyed by the gtest identifier that locks it."""
    values: dict[str, ExpectedValue] = {}

    book_model = OrnsteinUhlenbeck(mu=0.3, alpha=8.0, sigma=0.3)
    book_costs = TradingCosts(rate=0.05, cost=0.02)
    exponential_model = OrnsteinUhlenbeck(mu=1.3499, alpha=5.0, sigma=0.15)
    stop_loss_model = OrnsteinUhlenbeck(mu=0.5388, alpha=16.6677, sigma=0.1599)
    stop_loss_costs = TradingCosts(rate=0.05, cost=0.05)
    book_stop_loss = 0.04

    density_model = OrnsteinUhlenbeck(mu=0.998, alpha=0.0045, sigma=0.0038)
    density = hitting_time_density(1.02, model=density_model, first=1.04, second=1.0)
    values["HittingTimeDensityTest.hittingTimeDensityOutputTest"] = density
    values["OuModelTest.hittingTimeDensityOutputTest"] = density

    exit_level = solve_exit_level(
        lambda level: exit_residual(level, model=book_model, costs=book_costs),
        model=book_model,
        costs=book_costs,
    )
    values["OptimalTradingLevelsTest.optimalExitLevelOutputTest"] = exit_level
    values["TradingLevelsTest.exitLevelOutputTest"] = exit_level
    exit_level_exponential = solve_exit_level(
        lambda level: exit_residual_exponential(level, model=exponential_model, costs=book_costs),
        model=exponential_model,
        costs=book_costs,
    )
    values["OptimalTradingLevelsTest.optimalExitLevelExponentialOutputTest"] = (
        exit_level_exponential
    )
    values["TradingLevelsTest.exitLevelExponentialOutputTest"] = exit_level_exponential

    exit_level_stop_loss = solve_exit_level(
        lambda level: exit_residual_stop_loss(
            level, model=book_model, costs=book_costs, stop_loss=book_stop_loss
        ),
        model=book_model,
        costs=book_costs,
    )
    values["OptimalTradingLevelsTest.optimalExitLevelStopLossOutputTest"] = exit_level_stop_loss
    values["TradingLevelsTest.exitLevelStopLossOutputTest"] = exit_level_stop_loss

    entry_level_stop_loss = solve_level(
        lambda level: entry_residual_stop_loss(
            level,
            model=book_model,
            costs=book_costs,
            exit_level=0.455191,
            stop_loss=book_stop_loss,
        ),
        lower=book_stop_loss,
        upper=0.455191,
    )
    values["OptimalTradingLevelsTest.optimalEntryLevelStopLossOutputTest"] = entry_level_stop_loss
    values["TradingLevelsTest.entryLevelStopLossOutputTest"] = entry_level_stop_loss

    lower_entry_level_stop_loss = solve_level(
        lambda level: lower_entry_residual_stop_loss(
            level,
            model=book_model,
            costs=book_costs,
            exit_level=0.455191,
            stop_loss=book_stop_loss,
        ),
        lower=book_stop_loss,
        upper=0.13093,
    )
    values["OptimalTradingLevelsTest.optimalEntryLowerStopLossOutputTest"] = (
        lower_entry_level_stop_loss
    )
    values["TradingLevelsTest.entryLevelLowerBoundStopLossOutputTest"] = lower_entry_level_stop_loss

    entry_level = solve_level(
        lambda level: entry_residual(
            level, model=book_model, costs=book_costs, exit_level=0.466836
        ),
        lower=book_model.mu - 4.0 * book_model.stationary_deviation,
        upper=0.466836,
    )
    values["OptimalTradingLevelsTest.optimalEntryLevelOutputTest"] = entry_level
    values["TradingLevelsTest.entryLevelOutputTest"] = entry_level

    entry_level_exponential = solve_level(
        lambda level: entry_residual_exponential(
            level, model=exponential_model, costs=book_costs, exit_level=1.4093
        ),
        lower=exponential_model.mu - 4.0 * exponential_model.stationary_deviation,
        upper=1.4093,
    )
    values["OptimalTradingLevelsTest.optimalEntryLevelExponentialOutputTest"] = (
        entry_level_exponential
    )
    values["TradingLevelsTest.entryLevelExponentialOutputTest"] = entry_level_exponential

    chain_entry_level = solve_level(
        lambda level: entry_residual_stop_loss(
            level,
            model=book_model,
            costs=book_costs,
            exit_level=exit_level_stop_loss,
            stop_loss=book_stop_loss,
        ),
        lower=book_stop_loss,
        upper=exit_level_stop_loss,
    )
    chain_lower_entry_level = solve_level(
        lambda level: lower_entry_residual_stop_loss(
            level,
            model=book_model,
            costs=book_costs,
            exit_level=exit_level_stop_loss,
            stop_loss=book_stop_loss,
        ),
        lower=book_stop_loss,
        upper=chain_entry_level,
    )
    values["TradingLevelsTest.levelOrderingChainTest.b_star"] = exit_level
    values["TradingLevelsTest.levelOrderingChainTest.b_star_stop_loss"] = exit_level_stop_loss
    values["TradingLevelsTest.levelOrderingChainTest.d_star"] = chain_entry_level
    values["TradingLevelsTest.levelOrderingChainTest.a_star"] = chain_lower_entry_level

    narrow_model = OrnsteinUhlenbeck(mu=0.995, alpha=0.02, sigma=0.003)
    values["OptimalMeanReversionTest.methodFOutputTest"] = function_f(
        1.01, model=narrow_model, rate=0.05
    )
    values["OptimalMeanReversionTest.methodGOutputTest"] = function_g(
        0.2, model=book_model, rate=0.05
    )
    values["OptimalMeanReversionTest.methodBOutputTest"] = exit_residual(
        0.4, model=book_model, costs=book_costs
    )
    values["OptimalMeanReversionTest.methodBStopLossOutputTest"] = exit_residual_stop_loss(
        0.28, model=book_model, costs=book_costs, stop_loss=0.2
    )
    values["OptimalMeanReversionTest.methodDOutputTest"] = entry_residual(
        0.4, model=book_model, costs=book_costs, exit_level=0.46683583
    )
    values["OptimalMeanReversionTest.methodDStopLossOutputTest"] = entry_residual_stop_loss(
        0.4, model=book_model, costs=book_costs, exit_level=0.466836, stop_loss=0.1
    )
    values["OptimalMeanReversionTest.methodAStopLossOutputTest"] = lower_entry_residual_stop_loss(
        0.4, model=book_model, costs=book_costs, exit_level=0.466836, stop_loss=-0.3
    )
    values["OptimalMeanReversionTest.methodAboveVOutputTest"] = value_function(
        0.55, model=book_model, costs=book_costs, exit_level=0.466836
    )
    values["OptimalMeanReversionTest.methodBelowVOutputTest"] = value_function(
        0.15, model=book_model, costs=book_costs, exit_level=0.466836
    )
    values["OptimalMeanReversionTest.methodAboveVStopLossOutputTest"] = value_function_stop_loss(
        0.6,
        model=stop_loss_model,
        costs=stop_loss_costs,
        exit_level=0.567304,
        stop_loss=0.4834,
    )
    values["OptimalMeanReversionTest.methodBelowVStopLossOutputTest"] = value_function_stop_loss(
        0.5,
        model=stop_loss_model,
        costs=stop_loss_costs,
        exit_level=0.567304,
        stop_loss=0.4834,
    )
    boundary_offset = 1e-6
    values["OptimalMeanReversionTest.methodVStopLossBoundaryTest.at_exit"] = (
        value_function_stop_loss(
            0.567304 - boundary_offset,
            model=stop_loss_model,
            costs=stop_loss_costs,
            exit_level=0.567304,
            stop_loss=0.4834,
        )
    )
    values["OptimalMeanReversionTest.methodVStopLossBoundaryTest.at_stop_loss"] = (
        value_function_stop_loss(
            0.4834 + boundary_offset,
            model=stop_loss_model,
            costs=stop_loss_costs,
            exit_level=0.567304,
            stop_loss=0.4834,
        )
    )

    values["GaussianDistributionTest.cdfTest"] = float(norm.cdf(1.2, loc=0.996, scale=1.1))

    slow_model = OrnsteinUhlenbeck(mu=0.5, alpha=0.02, sigma=0.05)
    values["OrnsteinUhlenbeckModelTest.getUnconditionalVarianceOutputTest"] = (
        slow_model.stationary_deviation**2
    )
    values["OrnsteinUhlenbeckModelTest.getConditionalVarianceOutputTest.step_1"] = (
        slow_model.conditional_variance(1.0)
    )
    values["OrnsteinUhlenbeckModelTest.getConditionalVarianceOutputTest.step_2"] = (
        slow_model.conditional_variance(2.0)
    )
    zero_alpha_model = OrnsteinUhlenbeck(mu=0.5, alpha=0.0, sigma=0.05)
    near_zero_alpha_model = OrnsteinUhlenbeck(mu=0.5, alpha=1.1e-12, sigma=0.05)
    values[
        "OrnsteinUhlenbeckModelTest.getConditionalVarianceZeroAlphaLimitTest.zero_alpha_step_1"
    ] = zero_alpha_model.conditional_variance(1.0)
    values[
        "OrnsteinUhlenbeckModelTest.getConditionalVarianceZeroAlphaLimitTest.zero_alpha_step_2"
    ] = zero_alpha_model.conditional_variance(2.0)
    values[
        "OrnsteinUhlenbeckModelTest.getConditionalVarianceZeroAlphaLimitTest.near_zero_alpha_step_1"
    ] = near_zero_alpha_model.conditional_variance(1.0)

    values["HittingTimeOrnsteinUhlenbeckTest.hittingTimeCoreOutputTest"] = hitting_time_core(
        0.3, model=slow_model
    )
    values["HittingTimeOrnsteinUhlenbeckTest.optimalTradingFCoreOutputTest"] = float(
        0.1 ** (0.02 / slow_model.alpha - 1.0)
        * np.exp(slow_model.scaled_drift * (0.3 - slow_model.mu) * 0.1 - 0.1**2 / 2.0)
    )
    values["HittingTimeOrnsteinUhlenbeckTest.optimalTradingGCoreOutputTest"] = float(
        0.1 ** (0.02 / slow_model.alpha - 1.0)
        * np.exp(slow_model.scaled_drift * (slow_model.mu - 0.3) * 0.1 - 0.1**2 / 2.0)
    )
    values["HittingTimeOrnsteinUhlenbeckTest.optimalTradingLCoreOutputTest"] = level_l_star(
        model=book_model, costs=book_costs
    )

    linear_model = GeneralLinear(mu=-0.00143647, sigma=10.4573)
    values["GeneralLinearModelTest.GetUnconditionalVarianceTest"] = (
        linear_model.unconditional_variance
    )
    values["GeneralLinearModelTest.GetConditionalVarianceTest.step_1"] = (
        linear_model.conditional_variance(1.0)
    )
    values["GeneralLinearModelTest.GetConditionalVarianceTest.step_2"] = (
        linear_model.conditional_variance(2.0)
    )
    zero_mu_model = GeneralLinear(mu=0.0, sigma=0.05)
    near_zero_mu_model = GeneralLinear(mu=1.1e-12, sigma=0.05)
    values["GeneralLinearModelTest.GetConditionalVarianceZeroMuLimitTest.zero_mu_step_1"] = (
        zero_mu_model.conditional_variance(1.0)
    )
    values["GeneralLinearModelTest.GetConditionalVarianceZeroMuLimitTest.zero_mu_step_2"] = (
        zero_mu_model.conditional_variance(2.0)
    )
    values["GeneralLinearModelTest.GetConditionalVarianceZeroMuLimitTest.near_zero_mu_step_1"] = (
        near_zero_mu_model.conditional_variance(1.0)
    )

    entry_sums = ornstein_uhlenbeck_sums([0.5, 0.25, 0.5, 0.75, 1.5, 0.5])
    entry_estimate = ornstein_uhlenbeck_estimate(entry_sums)
    values["OuModelTest.ornsteinUhlenbeckMaximumLikelihoodOutputTest.mu"] = entry_estimate.mu
    values["OuModelTest.ornsteinUhlenbeckMaximumLikelihoodOutputTest.alpha"] = entry_estimate.alpha
    values["OuModelTest.ornsteinUhlenbeckMaximumLikelihoodOutputTest.sigma"] = entry_estimate.sigma

    fixed_sums = OrnsteinUhlenbeckSums(4.0, 3.5, 4.125, 3.375, 3.25, 6)
    fixed_estimate = ornstein_uhlenbeck_estimate(fixed_sums)
    values["OrnsteinUhlenbeckLikelihoodCalculateTest.ParameterTest.mu"] = fixed_estimate.mu
    values["OrnsteinUhlenbeckLikelihoodCalculateTest.ParameterTest.alpha"] = fixed_estimate.alpha
    values["OrnsteinUhlenbeckLikelihoodCalculateTest.ParameterTest.sigma"] = fixed_estimate.sigma
    updated_sums = update_sums(fixed_sums, new_observation=0.75, last_observation=1.0)
    updated_estimate = ornstein_uhlenbeck_estimate(updated_sums)
    values["OrnsteinUhlenbeckLikelihoodUpdateTest.ParameterTest.mu"] = updated_estimate.mu
    values["OrnsteinUhlenbeckLikelihoodUpdateTest.ParameterTest.alpha"] = updated_estimate.alpha
    values["OrnsteinUhlenbeckLikelihoodUpdateTest.ParameterTest.sigma"] = updated_estimate.sigma
    for field_name, field_value in vars(updated_sums).items():
        values[f"OrnsteinUhlenbeckLikelihoodUpdateTest.ComponentsTest.{field_name}"] = field_value
    series_sums = ornstein_uhlenbeck_sums([0.5, 0.25, 0.5, 0.75, 1.5, 1.0])
    for field_name, field_value in vars(series_sums).items():
        values[f"OrnsteinUhlenbeckLikelihoodCalculateTest.ComponentsTest.{field_name}"] = (
            field_value
        )

    values["SolverBoundsTest.upperSolverBoundOutputTest"] = (
        stop_loss_model.mu + 4.0 * stop_loss_model.stationary_deviation
    )
    values["SolverBoundsTest.lowerSolverBoundOutputTest"] = (
        stop_loss_model.mu - 4.0 * stop_loss_model.stationary_deviation
    )

    wide_model = OrnsteinUhlenbeck(mu=0.996, alpha=5.1, sigma=1.1)
    adaptive, _ = quad(
        lambda point: hitting_time_core(point, model=wide_model),
        0.8,
        1.05,
        epsabs=0.0,
        epsrel=1e-12,
    )
    values["AdaptiveIntegrationFunctionTest.OutputTest"] = adaptive
    values["SemiInfiniteIntegrationFunctionTest.OutputTest"] = kernel_integral(
        0.9, model=wide_model, rate=0.03, mean_sign=1.0, lower=0.8
    )
    values["AdaptiveCentralDifferentiationFunctionTest.OutputTest"] = 2.0
    values["BrentSolverFunctionTest.OutputTest"] = float(np.sqrt(5.0))

    values.update(filter_states_expected_values())
    values.update(kca_expected_values())
    return values


def snake_case(name: str) -> str:
    """Snake case of a camel case or pascal case gtest identifier."""
    spaced = re.sub(r"(?<=[a-z0-9])([A-Z])", r"_\1", name)
    spaced = re.sub(r"(?<=[A-Z])([A-Z])(?=[a-z])", r"_\1", spaced)
    return spaced.lower()


def format_value(value: ExpectedValue) -> str:
    """Listing representation of 1 expected value."""
    if isinstance(value, (list, int)):
        return repr(value)
    return repr(float(value))


def format_constant(name: str, value: ExpectedValue) -> str:
    """Inline constexpr definition of 1 expected value."""
    if isinstance(value, bool) or (isinstance(value, int) and value < 0):
        raise ValueError(f"{name} must be a count, not {value!r}.")
    if isinstance(value, list):
        if value and isinstance(value[0], list):
            if any(len(row) != len(value[0]) for row in value):
                raise ValueError(f"{name} has ragged rows.")
            rows = ", ".join(
                "{" + ", ".join(repr(float(element)) for element in row) + "}" for row in value
            )
            return (
                f"inline constexpr std::array<std::array<double, {len(value[0])}>, "
                f"{len(value)}> {name} = {{{{{rows}}}}};"
            )
        elements = ", ".join(repr(float(element)) for element in value)
        return f"inline constexpr std::array<double, {len(value)}> {name} = {{{elements}}};"
    if isinstance(value, int):
        return f"inline constexpr unsigned {name} = {value}u;"
    return f"inline constexpr double {name} = {repr(float(value))};"


def generate_header(values: dict[str, ExpectedValue]) -> str:
    """C++ header of the expected values, 1 nested namespace per test suite."""
    suites: dict[str, dict[str, ExpectedValue]] = {}
    for key, value in values.items():
        suite, _, remainder = key.partition(".")
        test_name, _, field = remainder.partition(".")
        constant = snake_case(test_name) + (f"_{field}" if field else "")
        namespace = suites.setdefault(snake_case(suite), {})
        if constant in namespace:
            raise ValueError(f"Duplicate constant {constant} in namespace {snake_case(suite)}.")
        namespace[constant] = value

    lines = [
        "// Generated by tools/derivations/expected_values.py. Do not edit.",
        "// clang-format off",
        "#ifndef STOCHASTIC_MODELS_TESTS_SUPPORT_EXPECTED_VALUES_H",
        "#define STOCHASTIC_MODELS_TESTS_SUPPORT_EXPECTED_VALUES_H",
        "",
        "#include <array>",
        "",
        "namespace expected {",
        "",
    ]
    for suite in sorted(suites):
        lines.append(f"namespace {suite} {{")
        for constant in sorted(suites[suite]):
            lines.append(format_constant(constant, suites[suite][constant]))
        lines.append("}")
        lines.append("")
    lines.append("} // namespace expected")
    lines.append("")
    lines.append("#endif // STOCHASTIC_MODELS_TESTS_SUPPORT_EXPECTED_VALUES_H")
    return "\n".join(lines) + "\n"


CHECK_RELATIVE_TOLERANCE = 1e-10
VALUE_LITERAL = re.compile(r"(?<=[ {])-?\d+(?:\.\d+)?(?:e[-+]?\d+)?(?=[,;}])")


def check_header(path: Path, values: dict[str, ExpectedValue]) -> list[str]:
    """Differences between the header at `path` and the generated header.

    A line may differ only in its numeric literals, and each literal only
    within the relative tolerance. A zero must stay zero.
    """
    committed = path.read_text(encoding="utf-8").splitlines()
    generated = generate_header(values).splitlines()
    failures: list[str] = []
    if len(committed) != len(generated):
        failures.append(f"{path}: {len(committed)} lines, generated {len(generated)}")
    for number, (old_line, new_line) in enumerate(zip(committed, generated), start=1):
        if old_line == new_line:
            continue
        if VALUE_LITERAL.sub("?", old_line) != VALUE_LITERAL.sub("?", new_line):
            failures.append(f"{path}:{number}: text differs")
            failures.append(f"  committed: {old_line}")
            failures.append(f"  generated: {new_line}")
            continue
        for old_text, new_text in zip(
            VALUE_LITERAL.findall(old_line), VALUE_LITERAL.findall(new_line)
        ):
            old_value, new_value = float(old_text), float(new_text)
            if not math.isclose(
                old_value, new_value, rel_tol=CHECK_RELATIVE_TOLERANCE, abs_tol=0.0
            ):
                failures.append(
                    f"{path}:{number}: {old_text} is not within "
                    f"{CHECK_RELATIVE_TOLERANCE:g} relative of {new_text}"
                )
    return failures


def main() -> None:
    """Print every expected value, write them as a C++ header, or check a header."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    output = parser.add_mutually_exclusive_group()
    output.add_argument(
        "--header",
        metavar="PATH",
        help="write the expected values as a C++ header at PATH instead of printing them",
    )
    output.add_argument(
        "--check",
        metavar="PATH",
        help="compare the header at PATH with the expected values and exit 1 on a difference",
    )
    arguments = parser.parse_args()
    values = expected_values()
    if arguments.header is not None:
        Path(arguments.header).write_text(generate_header(values), encoding="utf-8")
        return
    if arguments.check is not None:
        failures = check_header(Path(arguments.check), values)
        for failure in failures:
            print(failure, file=sys.stderr)
        if failures:
            sys.exit(1)
        print(
            f"{arguments.check}: {len(values)} values within {CHECK_RELATIVE_TOLERANCE:g} relative"
        )
        return
    for name, value in values.items():
        print(f"{name} {format_value(value)}")


if __name__ == "__main__":
    main()
