#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
	#define M_PI		3.14159265358979323846
#endif

///////////////////////////////////////////////////////////////////
// Vec
///////////////////////////////////////////////////////////////////

template <typename T>
struct Vec3{
	T x;
	T y;
	T z;
	Vec3(): x(0), y(0), z(0) {};
	Vec3(T x, T y, T z) : x(x), y(y), z(z) {};
	Vec3<T> normalized();
};

template <typename T>
Vec3<T> Vec3<T>::normalized() {
    double d = sqrt(x*x + y*y + z*z);
    return (d != 0.) ? Vec3<T>(x/d, y/d, z/d) : Vec3<T>();
}

template <typename T>
Vec3<T> operator * (const Vec3<T> &a, const Vec3<T> &b){
	return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

template <typename T>
Vec3<T> operator + (const Vec3<T> &a, const Vec3<T> &b){
	return {a.x + b.x, a.y + b.y, a.z + b.z};
}

template <typename T>
Vec3<T>& operator += (Vec3<T> &a, const Vec3<T> &b){
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}

template <typename T>
Vec3<T>& operator -= (Vec3<T> &a, const Vec3<T> &b){
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}

template <typename T>
Vec3<T> operator - (const Vec3<T> &a, const Vec3<T> &b){
	return {a.x - b.x, a.y - b.y, a.z - b.z};
}

template <typename T>
Vec3<T> operator - (const Vec3<T> &a){
	return {- a.x, - a.y, - a.z};
}

template <typename T>
Vec3<T> operator * (const double &p, const Vec3<T> &a){
	return {p*a.x, p*a.y, p*a.z};
}

template <typename T>
Vec3<T> operator * (const Vec3<T> &a, const double &p){
	return {p*a.x, p*a.y, p*a.z};
}

template <typename T>
Vec3<T> operator / (const Vec3<T> &a, const double &p){
	return {a.x/p, a.y/p, a.z/p};
}

template <typename T>
double dot(const Vec3<T> &a, const Vec3<T> &b){
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

template <typename T>
double norm(const Vec3<T> &a){
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

template <typename T>
void to_console(Vec3<T> a){
    printf("{%2.4e, %2.4e, %2.4e}", a.x, a.y, a.z);
}

typedef Vec3<double> Vector;

void testVector() {
    Vector a = Vector(1, 2, 3);
    Vector b = Vector(3, 2, 1);
    Vector c;
    to_console(a); printf("\n");
    to_console(b); printf("\n");
    to_console(a + b); printf("\n");
    to_console(a - b); printf("\n");
    to_console(- a); printf("\n");
    to_console(2*a); printf("\n");
    to_console(a*2); printf("\n");
    to_console(a/2); printf("\n");
    to_console(a*b); printf("\n");
    to_console(a.normalized()); printf("\n");
    printf("%e\n", norm(a));
    printf("%e\n", dot(a, b));
}
