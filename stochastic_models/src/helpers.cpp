#include "numeric_utils/helpers.h"

void ValuesSquared(std::vector<double>& vec) {
    for (double& item : vec) {
        item *= item;
    }
}
