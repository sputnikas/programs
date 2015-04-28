#ifndef SDE_H
#define SDE_H

#include <stdio.h>
#include <math.h>
#include <string.h>

void right_part(double *y, double t, int n, double *a);
void sum_vector(int n, double *y1, double a1, double *y2, double a2, double *result);
void euler(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result);
void runge2a(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result);
void runge2b(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result);
void runge2c(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result);
void runge3(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result);
void runge4a(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result);
void runge4b(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result);
void sde(double *y, double *y0, double t0, double tmax, int n, int N, void rp(double*, double, int, double*), int method);
double norm(double *y, double *y_p, int n);
void massive_copy(double *y, double *y_more, int N, int n, int i);
void sde_we(double *y, double *y0, double t0, double tmax, int n, int N, void rp(double*, double, int, double*), int method, double epsilon = 1e-6);

#define EULER        1 // Метод Эйлера
#define RUNGESA      2 // Метод Рунге-Кутты 2-го порядка
#define RUNGESB      3 // Метод Рунге-Кутты 2-го порядка
#define RUNGESC      4 // Метод Рунге-Кутты 2-го порядка
#define RUNGET       5 // Метод Рунге-Кутты 3-го порядка
#define RUNGEFA      6 // Метод Рунге-Кутты 4-го порядка
#define RUNGEFB      7 // Метод Рунге-Кутты 4-го порядка
#define ADAMSS       8 //Метод Адамса 2-го порядка
#define ADAMST       9 //Метод Адамса 3-го порядка
#define ADAMSF      10 //Метод Адамса 4-го порядка
#define ADAMSSE     11 //Метод Адамса 2-го порядка с точностью на шаге
#define ADAMSTE     12 //Метод Адамса 3-го порядка с точностью на шаге
#define ADAMSFE     13 //Метод Адамса 4-го порядка с точностью на шаге


#endif // SDE_H
