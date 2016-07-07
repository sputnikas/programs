#include <cmath>
#include <string>
#include <sstream>

#include "Array.h"
#include "Function.h"

template <typename T>
Array<T>::Array()
{
    size = 0;
    data = new T[1];
}

template <typename T>
Array<T>::Array(T *data, const uint32_t & size) : data(data), size(size)
{
}

template <typename T>
Array<T>::Array(const T &a, const uint32_t & size) : size(size)
{
    data = new T[size];
    for (uint32_t i = 0; i<size; i++)
    {
        data[i] = a;
    }
}

template <typename T>
Array<T>::Array(const uint32_t & size) : size(size)
{
    data = new T[size];
    for (uint32_t i = 0; i<size; i++)
    {
        data[i] = 0;
    }
}

template <typename T>
Array<T>::Array(const Array &a)
{
    size = a.size;
    data = new T[size];
    for (uint32_t i = 0; i<size; i++)
    {
        data[i] = a.data[i];
    }
}

template <typename T>
Array<T>::~Array()
{
    delete [] data;
}

template <typename T>
T Array<T>::get(const uint32_t &i)
{
    return (i < size) ? data[i] : 0;
}

template <typename T>
void Array<T>::set(const uint32_t &i, const T &t)
{
    if (i < size)
    {
        data[i] = t;
    }
}

template <typename T>
Array<T> Array<T>::operator +(const Array<T> &a)
{
    Array<T> result = (size > a.size) ? (*this) : a;
    uint32_t smin = (size > a.size) ? a.size : size;
    for (uint32_t i = 0; i<smin; i++)
    {
        result.data[i] += (size > a.size) ? a.data[i] : data[i];
    }
    return result;
}

template <typename T>
Array<T> Array<T>::operator -()
{
    Array<T> result = Array<T>(*this);
    for (uint32_t i = 0; i<size; i++)
    {
        result.data[i] = - data[i];
    }
    return result;
}

template <typename T>
Array<T> Array<T>::operator -(const Array<T> &a)
{
    Array<T> result = (size > a.size) ? (*this) : a;
    uint32_t smin = (size > a.size) ? a.size : size;
    for (uint32_t i = 0; i<smin; i++)
    {
        result.data[i] -= (size > a.size) ? a.data[i] : -data[i];
    }
    return result;
}

template <typename T>
Array<T> Array<T>::operator *(const T &a)
{
    Array<T> result = Array<T>(*this);
    for (uint32_t i = 0; i<size; i++)
    {
        result.data[i] *= a;
    }
    return result;
}

template <typename T>
Array<T> Array<T>::operator /(const T &a)
{
    Array<T> result = Array<T>(*this);
    for (uint32_t i = 0; i<size; i++)
    {
        result.data[i] *= a;
    }
    return result;
}

template<typename T>
std::ostream & operator << (std::ostream &out, const Array<T> &a)
{
    out << "{ ";
    if (a.size > 0)
    {
        for (uint32_t i = 0; i<a.size - 1; i++)
            out << a.data[i] << ", ";
        out << a.data[a.size - 1];
    }
    out << " }";
    return out;
}

void test_array()
{
    typedef Array<double> arr;
    arr a;
    std::cout << a << std::endl;
    arr b = arr(5);
    std::cout << b << std::endl;
    arr c = arr(2., 6);
    std::cout << c << std::endl;
    arr d = c*2.;
    std::cout << d << std::endl;
    std::cout << b + c << std::endl;
    std::cout << b - c << std::endl;
    std::cout << c/2 << std::endl;
}
