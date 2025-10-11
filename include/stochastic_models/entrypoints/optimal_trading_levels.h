#ifndef OPTIMAL_TRADING_LEVELS_H
#define OPTIMAL_TRADING_LEVELS_H

/**
 * @brief Calculate the optimal trading level exit value when a stop loss is
 * provided.
 *
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param stop_loss The stop loss value.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level exit value.
 */
const double optimalExitLevel(
    const double mu,
    const double alpha,
    const double sigma,
    const double stop_loss,
    const double r,
    const double c
);
/**
 * @brief Calculate the optimal trading level exit value.
 *
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level exit value.
 */
const double optimalExitLevel(
    const double mu,
    const double alpha,
    const double sigma,
    const double r,
    const double c
);
/**
 * @brief Calculate the optimal trading level exit value for an exponential
 * stochatic model.
 *
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level exit value.
 */
const double optimalExitLevelExponential(
    const double mu,
    const double alpha,
    const double sigma,
    const double r,
    const double c
);
/**
 * @brief Calculate the lower bound of the optimal trading entry interval for an
 * exponential stochastic model.
 *
 * @param d_star The optimal trading level upper entry value.
 * @param b_star The optimal trading level exit value.
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level entry value.
 */
const double optimalEntryLevelLowerExponential(
    const double d_star,
    const double b_star,
    const double mu,
    const double alpha,
    const double sigma,
    const double r,
    const double c
);
/**
 * @brief Calculate the lower bound of the optimal trading entry interval when a
 * stop loss is provided.
 *
 * @param d_star The optimal trading level upper entry value.
 * @param b_star The optimal trading level exit value.
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param stop_loss The stop loss value.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level entry value.
 */
const double optimalEntryLevelLower(
    const double d_star,
    const double b_star,
    const double mu,
    const double alpha,
    const double sigma,
    const double stop_loss,
    const double r,
    const double c
);
/**
 * @brief Calculate the optimal trading level entry value for an exponential
 * stochatsic model.
 *
 * @param b_star The optimal trading level exit value.
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level entry value.
 */
const double optimalEntryLevelExponential(
    const double b_star,
    const double mu,
    const double alpha,
    const double sigma,
    const double r,
    const double c
);
/**
 * @brief Calculate the optimal trading level entry value when a stop loss is
 * provided.
 *
 * @param b_star The optimal trading level exit value.
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param stop_loss The stop loss value.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level entry value.
 */
const double optimalEntryLevel(
    const double b_star,
    const double mu,
    const double alpha,
    const double sigma,
    const double stop_loss,
    const double r,
    const double c
);
/**
 * @brief Calculate the optimal trading level entry value.
 *
 * @param b_star The optimal trading level exit value.
 * @param mu The series mean.
 * @param alpha The series mean reversion speed.
 * @param sigma The series volatility.
 * @param r The risk free rate.
 * @param c The transaction cost.
 * @return const double The optimal trading level entry value.
 */
const double optimalEntryLevel(
    const double b_star,
    const double mu,
    const double alpha,
    const double sigma,
    const double r,
    const double c
);

#endif // OPTIMAL_TRADING_LEVELS_H
