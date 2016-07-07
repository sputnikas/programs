#include <iostream>

#include "Function.h"

template <typename T>
T min(T a, T b)
{
    return (a > b) ? b : a;
}

template <typename T>
T max(T a, T b)
{
    return (a > b) ? a : b;
}

void test_function()
{
    std::cout << "1 = " << min(1, 2) << std::endl;
    std::cout << "2 = " << max(1, 2) << std::endl;
}
