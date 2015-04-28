#ifndef PDE_H
#define PDE_H

#include <stdio.h>
#include <math.h>

void pdepe(int Method,
           double *u,
           double *up,
           int N,
           double *z,
           double t,
           double dt,
           double cdz2(double z, double t),
           double source(double z, double t),
           double c1,
           double dc1,
           double c2,
           double dc2,
           double b1,
           double b2
          );

void pdehe(int Method,
           double *u,
           double *up,
           double *upp,
           int N,
           double *z,
           double t,
           double dt,
           double cdz2(double z, double t),
           double source(double z, double t),
           double c1,
           double dc1,
           double c2,
           double dc2,
           double b1,
           double b2
          );

void go_throw(double *u,
              double *y,
              double dx,
              int k,
              int m,
              int ib,
              int M,
              int N,
              double x_c1,
              double x_dcx1,
              double x_dcy1,
              double x_b1(double y)
              );

void pdeee(double *u,
           int N,
           int M,
           double *x,
           double *y,
           double omega,
           double eps,
           double source(double x, double y, double u),
           double x_c1,
           double x_dcx1,
           double x_dcy1,
           double x_b1(double y),
           double x_c2,
           double x_dcx2,
           double x_dcy2,
           double x_b2(double y),
           double y_c1,
           double y_dcx1,
           double y_dcy1,
           double y_b1(double x),
           double y_c2,
           double y_dcx2,
           double y_dcy2,
           double y_b2(double x)
          );

void pdeee(double *u,
           int N,
           int M,
           double *x,
           double *y,
           double omega,
           double eps,
           int i1,
           int j1,
           double u1,
           int i2,
           int j2,
           double u2,
           double x_c1,
           double x_dcx1,
           double x_b1(double y),
           double x_c2,
           double x_dcx2,
           double x_b2(double y),
           double y_c1,
           double y_dcy1,
           double y_b1(double x),
           double y_c2,
           double y_dcy2,
           double y_b2(double x)
          );
#endif // PDE_H
