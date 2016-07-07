#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Vector3.h"

#ifndef M_PI
	#define M_PI		3.14159265358979323846
#endif

///////////////////////////////////////////////////////////////////
// Mat
///////////////////////////////////////////////////////////////////

template <typename T>
struct Mat3{
	T a[9];
	Mat3() {};
	Mat3(T a11, T a12, T a13, T a21, T a22, T a23, T a31, T a32, T a33)
    {
        a[0] = a11; a[1] = a12; a[2] = a13;
        a[3] = a21; a[4] = a22; a[5] = a23;
        a[6] = a31; a[7] = a32; a[8] = a33;
    };
    T& operator () (int i, int j) { return ((i<3) && (j<3)) ? a[i*3 + j] : a[0]; };
    const T& operator () (int i, int j) const { return ((i<3) && (j<3)) ? a[i*3 + j] : a[0]; };

    Mat3<T> transposed();
    Mat3<T> inversed();
    T A(int i, int j);
    T det();
};

template <typename T>
Mat3<T> Mat3<T>::transposed(){
    Mat3<T> r;
    for (int i = 0; i<3; i++)
        for (int j = 0; j<3; j++) {
            r(i, j) = a[j*3 + i];
        }
	return r;
}

template <typename T>
T Mat3<T>::A(int i, int j){
    int k1 = (i + 1) % 3, k2 = (i + 2) % 3;
    int n1 = (j + 1) % 3, n2 = (j + 2) % 3;
	return a[k1*3 + n1]*a[k2*3 + n2] - a[k1*3 + n2]*a[k2*3 + n1];
}

template <typename T>
T Mat3<T>::det(){
	return a[0]*A(0, 0) + a[1]*A(0, 1) + a[2]*A(0, 2);
}

template <typename T>
Mat3<T> Mat3<T>::inversed(){
    Mat3<T> r;
    T d = det();
    if (d != 0) {
        for (int i = 0; i<3; i++)
            for (int j = 0; j<3; j++)
                r(i, j) = A(j, i)/d;
    }
	return r;
}

template <typename T>
Mat3<T> operator * (const Mat3<T> &a, const Mat3<T> &b){
    Mat3<T> r;
    for (int i = 0; i<3; i++)
        for (int j = 0; j<3; j++) {
            r(i, j) = a(i, 0)*b(0, j) + a(i, 1)*b(1, j) + a(i, 2)*b(2, j);
        }
	return r;
}

template <typename T>
Vec3<T> operator * (const Mat3<T> &a, const Vec3<T> &b){
    Vec3<T> r;
    r.x = a(0, 0)*b.x + a(0, 1)*b.y + a(0, 2)*b.z;
    r.y = a(1, 0)*b.x + a(1, 1)*b.y + a(1, 2)*b.z;
    r.z = a(2, 0)*b.x + a(2, 1)*b.y + a(2, 2)*b.z;
	return r;
}

template <typename T>
Mat3<T> operator + (const Mat3<T> &a, const Mat3<T> &b){
    Mat3<T> r;
    for (int i = 0; i<3; i++)
        for (int j = 0; j<3; j++) {
            r(i, j) = a(i, j) + b(i, j);
        }
	return r;
}

template <typename T>
Mat3<T>& operator += (Mat3<T> &a, const Mat3<T> &b){
	return a + b;
}

template <typename T>
Mat3<T> operator - (const Mat3<T> &a, const Mat3<T> &b){
	Mat3<T> r;
    for (int i = 0; i<3; i++)
        for (int j = 0; j<3; j++) {
            r(i, j) = a(i, j) - b(i, j);
        }
	return r;
}

template <typename T>
Mat3<T>& operator -= (Mat3<T> &a, const Mat3<T> &b){
	return a - b;
}

template <typename T>
Mat3<T> operator - (const Mat3<T> &a){
    Mat3<T> r;
    for (int i = 0; i<3; i++)
        for (int j = 0; j<3; j++) {
            r(i, j) = - a(i, j);
        }
	return r;
}

template <typename T>
Mat3<T> operator * (const double &p, const Mat3<T> &a){
	Mat3<T> r;
    for (int i = 0; i<3; i++)
        for (int j = 0; j<3; j++) {
            r(i, j) = p*a(i, j);
        }
	return r;
}

template <typename T>
Mat3<T> operator * (const Mat3<T> &a, const double &p){
	Mat3<T> r;
    for (int i = 0; i<3; i++)
        for (int j = 0; j<3; j++) {
            r(i, j) = p*a(i, j);
        }
	return r;
}

template <typename T>
Mat3<T> operator / (const Mat3<T> &a, const double &p){
	Mat3<T> r;
    for (int i = 0; i<3; i++)
        for (int j = 0; j<3; j++) {
            r(i, j) = a(i, j) / p;
        }
	return r;
}

template <typename T>
void to_console(Mat3<T> a) {
    printf("{{%2.4e, %2.4e, %2.4e}, \n", a(0, 0), a(0, 1), a(0, 2));
    printf("{%2.4e, %2.4e, %2.4e}, \n",  a(1, 0), a(1, 1), a(1, 2));
    printf("{%2.4e, %2.4e, %2.4e}}, \n", a(2, 0), a(2, 1), a(2, 2));
}

typedef Mat3<double> Matrix;

void testMatrix() {
    Matrix a = Matrix(1, 2, 3, 4, 4, 4, 7, 8, 8);
    Matrix b = Matrix(1, 0, 0, 0, 2, 1, 0, 0, 3);
    Vector c = Vector(1, 0, 0);
    to_console(a);      printf("\n");
    to_console(b);      printf("\n");
    to_console(a + b);  printf("\n");
    to_console(a - b);  printf("\n");
    to_console(- a);    printf("\n");
    to_console(a * b);  printf("\n");
    to_console(a * c);  printf("\n");
    to_console(a / 2);  printf("\n");
    to_console(a * 2);  printf("\n");
    printf("%e\n", b.det());
    printf("%e\n", a.A(0, 1));
    to_console(a.inversed());  printf("\n");
    to_console(a.transposed());  printf("\n");
    to_console(a*a.inversed());  printf("\n");
}
