#include <iostream>
#include <chrono>
#include <vector>

template <typename T>
struct Complex {
    T x;
    T y;
    Complex() : x(0), y(0) {};
    Complex(T x, T y) : x(x), y(y) {};
    Complex(const T &x) : x(x), y(0) {};
};

typedef Complex<int>    CInt;
typedef Complex<float>  CFloat;
typedef Complex<double> CDouble;
typedef Complex<long double> CLDouble;


int main() {
    using namespace std::chrono;
    steady_clock::time_point t1, t2;

    int N = 1000000;
    //int *a = new int[N];
    //CInt *a = new CInt[N];
    //CFloat *a = new CFloat[N];
    //CDouble *a = new CDouble[N];
    double *a = new double[N];

    t1 = steady_clock::now();

    for (int i = 0; i<N; i++) {
        a[i] = i;
    }

    t2 = steady_clock::now();
    std::cout << "for a[i] \t";
    std::cout << (duration_cast<duration<double>>(t2 - t1)).count() << std::endl;

    t1 = steady_clock::now();

    for (int i = 0; i<N; i++) {
        *(a + i) = i;
    }

    t2 = steady_clock::now();
    std::cout << "for *(a + i) \t";
    std::cout << (duration_cast<duration<double>>(t2 - t1)).count() << std::endl;

    t1 = steady_clock::now();

    for (int i = 0; i<N; i++) {
        a[0] = i;
    }

    t2 = steady_clock::now();
    std::cout << "for a[0] \t";
    std::cout << (duration_cast<duration<double>>(t2 - t1)).count() << std::endl;

    t1 = steady_clock::now();

    for (int i = 0; i<N; i++) {
        a[N/2 + 1] = i;
    }

    t2 = steady_clock::now();
    std::cout << "for a[N - 1] \t";
    std::cout << (duration_cast<duration<double>>(t2 - t1)).count() << std::endl;
    return 0;
}
