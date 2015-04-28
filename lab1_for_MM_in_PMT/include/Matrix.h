#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

#ifndef NULL_MASSIVE
#define NULL_MASSIVE 0
#endif

//typedef float f_type;
class Matrix;
class LUP;

class Matrix{
public:
    int n_row;
    int n_column;
    float *a;
    bool lupped;
    int *p;

    Matrix();
    Matrix(int n_row, int n_column, float *elements, float scale);
    Matrix(int n_row, int n_column);
    ~Matrix();

    const Matrix &operator=(const Matrix &m);
    Matrix operator-();
    Matrix operator*(float scale);
    Matrix operator/(float scale);
    Matrix operator*(const Matrix &m);
    Matrix operator+(const Matrix &m);
    Matrix operator-(const Matrix &m);
    bool operator==(const Matrix &m);
    bool IsMatrix()const;
    bool IsQuad()const;
    void Transpose();
    std::string ToString(const char *ch = "\n")const;
    void Identity(int n);
    void Diagonal(int n, float *elements);
    void SwappingOfRows(int i, int j);
    void SwappingOfColumns(int i, int j);
    Matrix LowerTriangular()const;
    Matrix UpperTriangular()const;
    float A(int i, int j);
    void ToLUP();
    void FromLUP();
    Matrix L();
    Matrix U();
    Matrix P();
    float Det();
    void Wandermond(int n, float *elements);
    void Test();
};

Matrix operator*(const float scale, const Matrix &m);
void TestMatrix();
void TestLUP();
void TestLUP(Matrix &m);

#endif // MATRIX_H
