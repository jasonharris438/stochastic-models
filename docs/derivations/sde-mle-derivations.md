# SDE discretizations and MLE derivations

Standard results for the exact discretization and maximum-likelihood
estimation of the Ornstein–Uhlenbeck and general linear SDEs used in this
library. The time step is $\mathrm{d}t = 1$ throughout the likelihood layer
(implicit in all estimators).

**Notation.** For an $N$-point series there are $n = N - 1$ transition pairs.
The lagged values are $x_0 \dots x_{N-2}$ and the leads are $x_1 \dots x_{N-1}$.
The component sums are

$$S_x = \sum_i x_i, \quad S_y = \sum_i x_{i+1}, \quad S_{xx} = \sum_i x_i^2, \quad S_{yy} = \sum_i x_{i+1}^2, \quad S_{xy} = \sum_i x_i\,x_{i+1}.$$

## Ornstein–Uhlenbeck process

$$\mathrm{d}X = \alpha(\mu - X)\,\mathrm{d}t + \sigma\,\mathrm{d}W.$$

Exact transition over a step of length $t$:

$$X_{s+t} \mid X_s \sim \mathcal{N}\!\left(\mu + (X_s - \mu)\,e^{-\alpha t},\ \ \sigma^2\,\frac{1 - e^{-2\alpha t}}{2\alpha}\right).$$

The $\alpha \to 0$ limit of the variance is $\sigma^2 t$; the stationary
variance is $\dfrac{\sigma^2}{2\alpha}$.

Maximum-likelihood estimators with unit steps over the $n = N - 1$ pairs:

$$\hat{\mu} = \frac{S_y S_{xx} - S_x S_{xy}}{n\,(S_{xx} - S_{xy}) - \left(S_x^2 - S_x S_y\right)},$$

$$\hat{\alpha} = \ln\!\left(S_{xx} - 2\hat{\mu} S_x + n\hat{\mu}^2\right) - \ln\!\left(S_{xy} - \hat{\mu}(S_x + S_y) + n\hat{\mu}^2\right),$$

$$a = e^{-\hat{\alpha}},$$

$$\tilde{v} = \frac{1}{n}\left[\,S_{yy} - 2a\,S_{xy} + a^2 S_{xx} - 2\hat{\mu}(1 - a)(S_y - a S_x) + n\hat{\mu}^2 (1 - a)^2\,\right],$$

$$\hat{\sigma} = \sqrt{\tilde{v}\cdot \frac{2\hat{\alpha}}{1 - a^2}}.$$

Here $\tilde{v}$ estimates the transition variance $\sigma^2\,\dfrac{1 - a^2}{2\alpha}$;
the factor $\dfrac{2\alpha}{1 - a^2}$ converts it to the SDE diffusion
coefficient and has limit $1$ as $\alpha \to 0$. Note that $n$ is the **pair**
count, and the factor uses $2\alpha$, not $2e^{-\alpha}$.

## General linear SDE

$$\mathrm{d}X = \mu X\,\mathrm{d}t + \sigma\,\mathrm{d}W.$$

Exact transition over a step of length $t$:

$$X_{s+t} \mid X_s \sim \mathcal{N}\!\left(X_s\,e^{\mu t},\ \ \sigma^2\,\frac{e^{2\mu t} - 1}{2\mu}\right).$$

The $\mu \to 0$ limit of the variance is $\sigma^2 t$; a stationary variance
exists only for $\mu < 0$, equal to $\dfrac{\sigma^2}{2\lvert\mu\rvert}$.

Maximum-likelihood estimators with unit steps over the $n = N - 1$ pairs,
with $S_{xx}$ the **lagged** squares $x_0^2 \dots x_{N-2}^2$:

$$\hat{a} = \frac{S_{xy}}{S_{xx}}, \qquad \hat{\mu} = \ln \hat{a},$$

$$\tilde{v} = \frac{1}{n}\sum_i \left(x_{i+1} - \hat{a}\,x_i\right)^2,$$

$$\hat{\sigma} = \sqrt{\tilde{v}\cdot \frac{2\hat{\mu}}{e^{2\hat{\mu}} - 1}}.$$

The conversion factor has limit $1$ as $\mu \to 0$.

## Sampling at spacing $t$ (round-trip identities used by the tests)

Applying the unit-step estimator to data sampled every $t$ time units of the
continuous process recovers, exactly:

- **OU:** $\mu$ unchanged, $\hat{\alpha} = \alpha t$, $\hat{\sigma} = \sigma\sqrt{t}$.
- **GL:** $\hat{\mu} = \mu t$, $\hat{\sigma} = \sigma\sqrt{t}$.
