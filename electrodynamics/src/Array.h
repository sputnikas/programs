#pragma once

template <typename T>
class Array {
    unsigned int n;
    T *data;
public:
    Array();
    Array(unsigned int n);
    Array(T *a, unsigned int n);
    Array(unsigned int n, T a, ...);
    Array(const Array<T> &a);

    void sort();
    void add(const T &a);
    void remove(const unsigned int &k);

};

template <typename T> Array<T> operator + (Array<T> a, Array<T> b);
template <typename T> Array<T> operator - (Array<T> a, Array<T> b);
template <typename T> Array<T> operator - (Array<T> a);
template <typename T> Array<T> operator * (T a, Array<T> b);
template <typename T> Array<T> operator * (Array<T> a, T b);
template <typename T> Array<T> operator / (Array<T> a, T b);

template <typename T> Array<T> cart(Array<T> a, Array<T> b);
