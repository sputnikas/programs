#pragma once
#include "Types.h"
#include <cmath>
#include <cstdarg>
#include <iostream>

template <typename T, u32 N>
class Vector {
private:
    T data[N] = {0};
public:
    Vector(){};
    Vector(const T (&a)[N]);
    Vector(const T a, ...);

    Vector operator + (const Vector& a) const;
    Vector operator - (const Vector& a) const;
    Vector operator * (const T& a) const;
    Vector operator / (const T& a) const;
    Vector operator += (const Vector& a);
    Vector operator -= (const Vector& a);
    Vector operator *= (const T& a);
    Vector operator /= (const T& a);
    Vector operator - ();
    void toConsole(const char* ch = "\n");
};

template <typename T, u32 N> Vector<T, N>::Vector(const T (&a)[N]) {
    for (u32 i = 0; i<N; i++) {
        data[i] = a[i];
    }
}

template <typename T, u32 N>Vector<T, N>::Vector(const T a, ...) {
    va_list args;
    va_start(args, a);
    data[0] = a;
    for (u32 i = 1; i<N; i++) {
        data[i] = va_arg(args, T);
    }
    va_end(args);
}

template <typename T, u32 N>
Vector<T, N> Vector<T, N>::operator + (const Vector<T, N>& a) const {
    Vector<T, N> r;
    for (u32 i = 0; i<N; i++) {
        r.data[i] = data[i] + a.data[i];
    }
    return r;
}

template <typename T, u32 N>
void Vector<T, N>::toConsole(const char* ch) {
    std::cout << "{";
    for (u32 i = 0; i<N - 1; i++) {
        std::cout << data[i] << ", ";
    }
    std::cout << data[N - 1] << "}" << ch;
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

void testVector() {
    Vector<double, 3> a;
    a.toConsole();
    Vector<double, 3> b = Vector<double, 3>(1.0, 2.0, 3.0);
    b.toConsole();
    a = Vector<double, 3>({3.0, 2.0, 1.0});
    a.toConsole();
    Vector<double, 3> c = a + b;
    a.toConsole();
    b.toConsole();
    c.toConsole();

}
