#ifndef ARRAY_H
#define ARRAY_H
#include <cstdint>
#include <iostream>

template<typename T> class Array
{
public:
    T *data;
    uint32_t size;

    Array();
    Array(T *data, const uint32_t & size);
    Array(const T &a, const uint32_t & size);
    Array(const uint32_t & size);
    Array(const Array &a);
    ~Array();

    T get(const uint32_t &i);
    void set(const uint32_t &i, const T &t);

    Array operator -();
    Array operator +(const Array &a);
    Array operator -(const Array &a);
    Array operator *(const T &a);
    Array operator /(const T &a);

};

template<typename T>
std::ostream & operator << (std::ostream &out, const Array<T> &a);

void test_array();

#endif // ARRAY_H
