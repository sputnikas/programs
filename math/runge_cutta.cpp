#include <stdio.h>
#include <math.h>
#include <stdarg.h>
#include <iostream>

struct vec {
    double *data;
    int n;

    vec();
    vec(const int &n);
    vec(double *a, const int &n);
    vec(double a, const int &n);
    vec(const vec &a);
    vec(double *a, const double &b, const int &n);
    vec(double *a, double *b, const double &ka, const double &kb, const int &n);
    ~vec();
};

vec::vec() {
    n = 0;
    data = new double[1];
}

vec::vec(const int &n) : n(n){
    data = new double[n];
}

vec::vec(double *a, const int &n) : n(n){
    data = new double[n];
    for (int i = 0; i<n; i++){
        data[i] = a[i];
    }
}

vec::vec(double a, const int &n) : n(n){
    data = new double[n];
    for (int i = 0; i<n; i++){
        data[i] = a;
    }
}

vec::vec(const vec &a){
    n = a.n;
    data = new double[n];
    for (int i = 0; i<n; i++){
        data[i] = a.data[i];
    }
}

vec::vec(double *a, const double &b, const int &n) : n(n){
    data = new double[n];
    for (int i = 0; i<n; i++){
        data[i] = b*a[i];
    }
}

vec::vec(double *a, double *b, const double &ka, const double &kb, const int &n) : n(n){
    data = new double[n];
    for (int i = 0; i<n; i++){
        data[i] = ka*a[i] + kb*b[i];
    }
}

vec::~vec(){
    delete [] data;
}

vec operator + (const vec &a, const vec &b){
    return (a.n == b.n) ? vec(a.data, b.data, 1., 1., a.n) : vec();
}

vec operator - (const vec &a, const vec &b){
    return (a.n == b.n) ? vec(a.data, b.data, 1., -1., a.n) : vec();
}

vec operator * (const vec &a, const double &b){
    return vec(a.data, b, a.n);
}

vec operator * (const double &b, const vec &a){
    return vec(a.data, b, a.n);
}

vec operator / (const vec &a, const double &b){
    return vec(a.data, 1./b, a.n);
}

double norm(const vec &a){
    double s = 0.;
    for (int i = 0; i<a.n; i++){
        s += a.data[i]*a.data[i];
    }
    return sqrt(s);
}

std::ostream & operator << (std::ostream &out, const vec &a)
{
    out << "{ ";
    if (a.n > 0)
    {
        for (int i = 0; i<a.n - 1; i++)
            out << a.data[i] << ", ";
        out << a.data[a.n - 1];
    }
    out << " }";
    return out;
}

vec runge(vec y, double x, double h, vec (*f)(double, vec)){
    vec k1 = f(x, y);
    vec k2 = f(x + h/2, y + k1*h/2);
    vec k3 = f(x + h/2, y + k2*h/2);
    vec k4 = f(x + h, y + k3*h);
    return y + (k1 + 2.*k2 + 2.*k3 + k4)*h / 6;
}

double a, b;

vec f(double x, vec y){
    vec r = vec(2);
    r.data[0] = y.data[1];
    r.data[1] = - (a - 2*b*cos(2*x))*y.data[0];
    return r;
}

void test_vec(){
    vec a;
    vec b = vec(2);
    vec c = vec(2, 2);
    vec d = 2*c;
    double D[2] = {1, 2};
    vec f = vec(D, 2);
    vec g = vec(D, 2, 2);
    vec h = vec(D, D, 0.5, 0.1, 2);

    std::cout << a << std::endl;
    std::cout << b << std::endl;
    std::cout << c << std::endl;
    std::cout << d << std::endl;
    std::cout << f << std::endl;
    std::cout << g << std::endl;
    std::cout << h << std::endl;
}

int main() {
    a = 0.5;
    b = 0.25;
    double h = 0.005;
    int N = 8000;

    FILE *fp;
    fp = fopen("2.gp", "w");
    fprintf(fp, "plot '-' using 1:2 w l\n");

    double init[2] = {1., 0.};
    vec y0 = vec(2), y = vec(init, 2);
    for (int i = 0; i<N; i++){
        y0 = y;
        y = runge(y0, i*h, h, f);
        fprintf(fp, "%e %e %e\n", i*h, y.data[0], y.data[1]);
    }

    fprintf(fp, "\n");

    double init2[2] = {0., 1.};
    y = vec(init2, 2);
    for (int i = 0; i<N; i++){
        y0 = y;
        y = runge(y0, i*h, h, f);
        fprintf(fp, "%e %e %e\n", i*h, y.data[0], y.data[1]);
    }
    fclose(fp);

    return 0;
}
