#pragma once

#include "Types.h"
#include <iostream>
#include <cmath>

template <typename T> T factorial(const T &a, const T &b) {
    return (a > 0) ? factorial<T>(a - 1, b*a) : b;
}

template <typename T> T factorial(const T &a) {
    return factorial<T>(a, 1);
}

template <typename T> T factorial2(const T &a) {
    return (a > 0) ? a*factorial2<T>(a - 1) : 1;
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

#include <chrono>

void testMathFunction(){
    std::cout << "0! = " << factorial(0) <<  std::endl;
    std::cout << "1! = " << factorial(1) <<  std::endl;
    std::cout << "2! = " << factorial(2) <<  std::endl;
    std::cout << "4! = " << factorial(4) <<  std::endl;
    std::cout << "10! = " << factorial(10) <<  std::endl;

    using namespace std::chrono;
    steady_clock::time_point t1, t2;
    u32 N = 10000;
    f64 result;

    result = 0;
    t1 = steady_clock::now();
    for (u32 i = 0; i<N; i++) {
        result += factorial(15.);
    }
    t2 = steady_clock::now();
    std::cout << "duration " << N << " times \t";
    std::cout << (duration_cast<duration<double>>(t2 - t1)).count() << std::endl;
    std::cout << "sum result = " << result <<  std::endl;

    result = 0;
    t1 = steady_clock::now();
    for (u32 i = 0; i<N; i++) {
        result += factorial2(15.);
    }
    t2 = steady_clock::now();
    std::cout << "duration " << N << " times \t";
    std::cout << (duration_cast<duration<double>>(t2 - t1)).count() << std::endl;
    std::cout << "sum result = " << result <<  std::endl;
}
