#include <iostream>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp; // using namespace for plt
int main(){
    plt::plot({1,3,1,4});
    plt::show(); // display
}