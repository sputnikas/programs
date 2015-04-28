#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <ctime>
#include <ratio>
#include <chrono>

#define EPSILON 1e-40
#define PI 3.1415926
#define LEFT -2
#define RIGHT -1
#define MIDDLE 1
#define SIMPSONE 2
#define THIRD 3

void SwapRows(int i, int j, double *matrix, int n);
void LUP(int n, double *matrix, double *l, double *u, int *p);
std::string MatrixToString(int n, double *matrix);
void TestLUP();
void SLAE(int n, double *matrix, double *right_part, double *solve);
template <class X> std::string MassiveToString(int n, double *matrix);
double func(double x);
double Integrate(int n, double *nodes, double *weights);
double IntegrateWithStep(int N, double a, double b, int n, double *weights);
double IntegrateLagranje(int n, double a, double b, double *x);
void TestIntegrate();
void TestLejandre(int n);

#endif // INTEGRATE_H
