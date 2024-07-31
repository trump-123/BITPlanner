#include "BITStarforssl.h"

int main() {
    std::pair<double, double> x_start = {18, 8};  // Starting node
    std::pair<double, double> x_goal = {37, 18};  // Goal node
    double eta = 2;
    int iter_max = 200;

    BITStar bit(x_start, x_goal, eta, iter_max);
    bit.plan();

    return 0;
}