"""Independent replica of the expected values locked by the unit tests.

Every value is computed from the model formulas with scipy, not from the
library. Quadrature uses QAWS for the algebraic endpoint singularity and
brentq for the roots. Run with:

    uv run --with numpy --with scipy tools/derivations/expected_values.py
"""

from dataclasses import dataclass
from typing import Callable, Dict, List

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.stats import norm


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


def hitting_time_core(x: float, *, model: OrnsteinUhlenbeck) -> float:
    """Unnormalised hitting time kernel exp(x alpha (x - 2 mu) / sigma^2)."""
    return float(np.exp(x * model.alpha * (x - 2.0 * model.mu) / model.sigma**2))


def hitting_time_density(
    x: float, *, model: OrnsteinUhlenbeck, first: float, second: float
) -> float:
    """Ratio of the kernel integral over [second, x] to that over [second, first].

    The exponent is shifted by its value at `second` so that both integrals
    stay in floating point range.
    """
    shift = x_exponent(second, model=model)

    def shifted_kernel(point: float) -> float:
        return float(np.exp(x_exponent(point, model=model) - shift))

    numerator, _ = quad(shifted_kernel, second, x, epsabs=0.0, epsrel=1e-12)
    denominator, _ = quad(shifted_kernel, second, first, epsabs=0.0, epsrel=1e-12)
    return numerator / denominator


def x_exponent(point: float, *, model: OrnsteinUhlenbeck) -> float:
    """Exponent of the hitting time kernel at one point."""
    return point * model.alpha * (point - 2.0 * model.mu) / model.sigma**2


def kernel_integral(
    x: float,
    *,
    model: OrnsteinUhlenbeck,
    rate: float,
    mean_sign: float,
    power_shift: int = 0,
    lower: float = 0.0,
) -> float:
    """Integral of u^(r/alpha - 1 + shift) exp(sign k (x - mu) u - u^2 / 2) over [lower, inf).

    With sign = +1 this is the Leung-Li F function, with sign = -1 the G
    function. A power shift of one gives the derivative integrand. The
    algebraic singularity at zero is handled by QAWS on [0, 1].
    """
    power = rate / model.alpha - 1.0 + power_shift
    drift = mean_sign * model.scaled_drift * (x - model.mu)

    def smooth_part(u: float) -> float:
        return float(np.exp(drift * u - 0.5 * u * u))

    def full_integrand(u: float) -> float:
        return u**power * smooth_part(u)

    if lower > 0.0:
        tail, _ = quad(full_integrand, lower, np.inf, epsabs=0.0, epsrel=1e-12, limit=200)
        return tail
    head, _ = quad(smooth_part, 0.0, 1.0, weight="alg", wvar=(power, 0.0), epsabs=0.0, epsrel=1e-12)
    tail, _ = quad(full_integrand, 1.0, np.inf, epsabs=0.0, epsrel=1e-12, limit=200)
    return head + tail


def function_f(x: float, *, model: OrnsteinUhlenbeck, rate: float) -> float:
    """Leung-Li F(x)."""
    return kernel_integral(x, model=model, rate=rate, mean_sign=1.0)


def function_g(x: float, *, model: OrnsteinUhlenbeck, rate: float) -> float:
    """Leung-Li G(x)."""
    return kernel_integral(x, model=model, rate=rate, mean_sign=-1.0)


def derivative_f(x: float, *, model: OrnsteinUhlenbeck, rate: float) -> float:
    """F'(x) by differentiation under the integral sign."""
    return model.scaled_drift * kernel_integral(
        x, model=model, rate=rate, mean_sign=1.0, power_shift=1
    )


def derivative_g(x: float, *, model: OrnsteinUhlenbeck, rate: float) -> float:
    """G'(x) by differentiation under the integral sign."""
    return -model.scaled_drift * kernel_integral(
        x, model=model, rate=rate, mean_sign=-1.0, power_shift=1
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
    x: float, *, model: OrnsteinUhlenbeck, costs: TradingCosts, exit_level: float
) -> float:
    """Value function (b - c) F(x) / F(b) below the exit, x - c above it."""
    if x >= exit_level:
        return x - costs.cost
    rate = costs.rate
    return (
        (exit_level - costs.cost)
        * function_f(x, model=model, rate=rate)
        / function_f(exit_level, model=model, rate=rate)
    )


def solve_exit_level(
    residual: Callable[[float], float], *, model: OrnsteinUhlenbeck, costs: TradingCosts
) -> float:
    """Root of the residual inside the library's solver bracket."""
    lower = max(level_l_star(model=model, costs=costs), costs.cost)
    upper = model.mu + 4.0 * model.stationary_deviation
    return float(brentq(residual, lower, upper, xtol=1e-14, rtol=1e-14, maxiter=500))


def ornstein_uhlenbeck_sums(series: List[float]) -> OrnsteinUhlenbeckSums:
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
    """Sufficient statistics after one more lag and lead pair."""
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


@dataclass(frozen=True)
class KcaState:
    """Kinetic components state: three-state constant-acceleration Kalman filter."""

    transition_matrix: np.ndarray
    transition_covariance: np.ndarray
    current_state_mean: np.ndarray
    current_state_covariance: np.ndarray
    observation_matrix: np.ndarray


def kca_initial_state(series: List[float], *, step: float, noise: float) -> KcaState:
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


def kca_update(state: KcaState, *, observation: float, innovation_sigma: float) -> KcaState:
    """One predict and correct step of the Kalman filter."""
    transition = state.transition_matrix
    predicted_mean = transition @ state.current_state_mean
    predicted_covariance = (
        transition @ state.current_state_covariance @ transition.T + state.transition_covariance
    )
    observation_matrix = state.observation_matrix
    innovation_covariance = (
        observation_matrix @ predicted_covariance @ observation_matrix.T + innovation_sigma**2
    )
    gain = predicted_covariance @ observation_matrix.T @ np.linalg.inv(innovation_covariance)
    innovation = observation - (observation_matrix @ predicted_mean)[0]
    corrected_mean = predicted_mean + gain[:, 0] * innovation
    corrected_covariance = (np.eye(3) - gain @ observation_matrix) @ predicted_covariance
    return KcaState(
        transition_matrix=transition,
        transition_covariance=state.transition_covariance,
        current_state_mean=corrected_mean,
        current_state_covariance=corrected_covariance,
        observation_matrix=observation_matrix,
    )


def kca_expected_values() -> Dict[str, float]:
    """Expected numeric fields of the two KCA entry point tests."""
    values: Dict[str, float] = {}
    series = [
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
    initial = kca_initial_state(series, step=1.0, noise=0.001)
    values["kca_filter_test.getInitializedKcaStateTest.current_state_mean_0"] = float(
        initial.current_state_mean[0]
    )
    for row in range(3):
        for column in range(3):
            values[
                f"kca_filter_test.getInitializedKcaStateTest.transition_covariance_{row}{column}"
            ] = float(initial.transition_covariance[row, column])
    prior = KcaState(
        transition_matrix=np.array(
            [[1.0011961162353782, 1.0, 0.5], [0.0, 1.0, 1.0], [0.0, 0.0, 1.0]]
        ),
        transition_covariance=np.diag([0.12695229227341848, 0.001, 0.001]),
        current_state_mean=np.array([10.288741828687053, 0.0, 0.0]),
        current_state_covariance=np.zeros((3, 3)),
        observation_matrix=np.array([[1.0, 0.0, 0.0]]),
    )
    posterior = kca_update(prior, observation=10.3, innovation_sigma=0.1)
    values["kca_filter_test.getUpdatedKcaStateTest.current_state_mean_0"] = float(
        posterior.current_state_mean[0]
    )
    for row in range(3):
        for column in range(3):
            values[
                f"kca_filter_test.getUpdatedKcaStateTest.current_state_covariance_{row}{column}"
            ] = float(posterior.current_state_covariance[row, column])
    return values


def expected_values() -> Dict[str, float]:
    """Every expected value, keyed by test file and test name."""
    values: Dict[str, float] = {}

    book_model = OrnsteinUhlenbeck(mu=0.3, alpha=8.0, sigma=0.3)
    book_costs = TradingCosts(rate=0.05, cost=0.02)
    exponential_model = OrnsteinUhlenbeck(mu=1.3499, alpha=5.0, sigma=0.15)
    stop_loss_model = OrnsteinUhlenbeck(mu=0.5388, alpha=16.6677, sigma=0.1599)

    density_model = OrnsteinUhlenbeck(mu=0.998, alpha=0.0045, sigma=0.0038)
    density = hitting_time_density(1.02, model=density_model, first=1.04, second=1.0)
    values["hitting_time_test.hittingTimeDensityOutputTest"] = density
    values["ou_model_test.hittingTimeDensityOutputTest"] = density

    exit_level = solve_exit_level(
        lambda level: exit_residual(level, model=book_model, costs=book_costs),
        model=book_model,
        costs=book_costs,
    )
    values["optimal_trading_levels_test.optimalExitLevelOutputTest"] = exit_level
    values["trading_levels_test.exitLevelOutputTest"] = exit_level
    exit_level_exponential = solve_exit_level(
        lambda level: exit_residual_exponential(level, model=exponential_model, costs=book_costs),
        model=exponential_model,
        costs=book_costs,
    )
    values["optimal_trading_levels_test.optimalExitLevelExponentialOutputTest"] = (
        exit_level_exponential
    )
    values["trading_levels_test.exitLevelExponentialOutputTest"] = exit_level_exponential

    narrow_model = OrnsteinUhlenbeck(mu=0.995, alpha=0.02, sigma=0.003)
    values["optimal_mean_reversion_test.methodFOutputTest"] = function_f(
        1.01, model=narrow_model, rate=0.05
    )
    values["optimal_mean_reversion_test.methodGOutputTest"] = function_g(
        0.2, model=book_model, rate=0.05
    )
    values["optimal_mean_reversion_test.methodBOutputTest"] = exit_residual(
        0.4, model=book_model, costs=book_costs
    )
    values["optimal_mean_reversion_test.methodBStopLossOutputTest"] = exit_residual_stop_loss(
        0.28, model=book_model, costs=book_costs, stop_loss=0.2
    )
    values["optimal_mean_reversion_test.methodAboveVOutputTest"] = value_function(
        0.55, model=book_model, costs=book_costs, exit_level=0.466836
    )
    values["optimal_mean_reversion_test.methodBelowVOutputTest"] = value_function(
        0.15, model=book_model, costs=book_costs, exit_level=0.466836
    )
    values["optimal_mean_reversion_test.methodAboveVStopLossOutputTest"] = value_function(
        0.6, model=stop_loss_model, costs=TradingCosts(rate=0.05, cost=0.05), exit_level=0.567304
    )

    values["gaussian_distribution_test.cdfTest"] = float(norm.cdf(1.2, loc=0.996, scale=1.1))

    slow_model = OrnsteinUhlenbeck(mu=0.5, alpha=0.02, sigma=0.05)
    values["ornstein_uhlenbeck_test.getUnconditionalVarianceOutputTest"] = (
        slow_model.stationary_deviation**2
    )
    values["ornstein_uhlenbeck_test.getMeanOutputTest"] = slow_model.mu
    values["ornstein_uhlenbeck_test.hittingTimeCoreOutputTest"] = hitting_time_core(
        0.3, model=slow_model
    )
    values["ornstein_uhlenbeck_test.optimalTradingFCoreOutputTest"] = float(
        0.1 ** (0.02 / slow_model.alpha - 1.0)
        * np.exp(slow_model.scaled_drift * (0.3 - slow_model.mu) * 0.1 - 0.1**2 / 2.0)
    )
    values["ornstein_uhlenbeck_test.optimalTradingGCoreOutputTest"] = float(
        0.1 ** (0.02 / slow_model.alpha - 1.0)
        * np.exp(slow_model.scaled_drift * (slow_model.mu - 0.3) * 0.1 - 0.1**2 / 2.0)
    )
    values["ornstein_uhlenbeck_test.optimalTradingLCoreOutputTest"] = level_l_star(
        model=book_model, costs=book_costs
    )

    entry_sums = ornstein_uhlenbeck_sums([0.5, 0.25, 0.5, 0.75, 1.5, 0.5])
    entry_estimate = ornstein_uhlenbeck_estimate(entry_sums)
    values["ou_model_test.ornsteinUhlenbeckMaximumLikelihoodOutputTest.mu"] = entry_estimate.mu
    values["ou_model_test.ornsteinUhlenbeckMaximumLikelihoodOutputTest.alpha"] = (
        entry_estimate.alpha
    )
    values["ou_model_test.ornsteinUhlenbeckMaximumLikelihoodOutputTest.sigma"] = (
        entry_estimate.sigma
    )

    fixed_sums = OrnsteinUhlenbeckSums(4.0, 3.5, 4.125, 3.375, 3.25, 6)
    fixed_estimate = ornstein_uhlenbeck_estimate(fixed_sums)
    values["ornstein_uhlenbeck_likelihood_test.CalculateParameterTest.mu"] = fixed_estimate.mu
    values["ornstein_uhlenbeck_likelihood_test.CalculateParameterTest.alpha"] = fixed_estimate.alpha
    values["ornstein_uhlenbeck_likelihood_test.CalculateParameterTest.sigma"] = fixed_estimate.sigma
    updated_sums = update_sums(fixed_sums, new_observation=0.75, last_observation=1.0)
    updated_estimate = ornstein_uhlenbeck_estimate(updated_sums)
    values["ornstein_uhlenbeck_likelihood_test.UpdateParameterTest.mu"] = updated_estimate.mu
    values["ornstein_uhlenbeck_likelihood_test.UpdateParameterTest.alpha"] = updated_estimate.alpha
    values["ornstein_uhlenbeck_likelihood_test.UpdateParameterTest.sigma"] = updated_estimate.sigma
    for field_name, field_value in vars(updated_sums).items():
        values[f"ornstein_uhlenbeck_likelihood_test.UpdateComponentsTest.{field_name}"] = (
            field_value
        )
    series_sums = ornstein_uhlenbeck_sums([0.5, 0.25, 0.5, 0.75, 1.5, 1.0])
    for field_name, field_value in vars(series_sums).items():
        values[f"ornstein_uhlenbeck_likelihood_test.CalculateComponentsTest.{field_name}"] = (
            field_value
        )

    values["utils_test.upperSolverBoundOutputTest"] = (
        stop_loss_model.mu + 4.0 * stop_loss_model.stationary_deviation
    )
    values["utils_test.lowerSolverBoundOutputTest"] = (
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
    values["utils_test.AdaptiveIntegrationOutputTest"] = adaptive
    values["utils_test.SemiInfiniteIntegrationOutputTest"] = kernel_integral(
        0.9, model=wide_model, rate=0.03, mean_sign=1.0, lower=0.8
    )
    values["utils_test.AdaptiveCentralDifferentiationOutputTest"] = 2.0
    values["utils_test.BrentSolverOutputTest"] = float(np.sqrt(5.0))
    values.update(kca_expected_values())
    return values


def main() -> None:
    """Print every expected value with full double precision."""
    for name, value in expected_values().items():
        print(f"{name} {float(value)!r}")


if __name__ == "__main__":
    main()
