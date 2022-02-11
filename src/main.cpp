#include "ornstein_uhlenbeck.h"

#include <iostream>

int main()
{
    OrnsteinUhlenbeckModel *model = new OrnsteinUhlenbeckModel();
    std::vector<const double> vec = (*model).Simulate(1000);

    for(size_t i{} ; i < vec.size() ; i++){
    std::cout << vec[i] << std::endl;
    }
    return 0;
}