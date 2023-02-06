#ifndef _LIKELIHOOD_H_
#define _LIKELIHOOD_H_
#include <map>
#include <string>
#include <vector>
/**
 * Abstract base likelihood calculator class.
 */
class LikelihoodCalculatuor {
   protected:
    /**
     * @brief Empty pointer to a map that stores likelihood calculation
     * components used in generating process parameters.
     *
     */
    const std::map<std::string, const double>* components;
    LikelihoodCalculatuor();

   public:
    /**
     * Calculates and returns the maximised likelihood parameters.
     *
     * @param data The data series to calculate the likelihood of.
     * @returns Const reference to mean of Gaussian Distribution.
     */
    virtual const std::map<std::string, const double> calculate(
        const std::vector<double>& data) = 0;
    ~LikelihoodCalculatuor();
};
#endif  // _LIKELIHOOD_H_
