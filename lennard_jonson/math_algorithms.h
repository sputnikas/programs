#ifndef MATH_ALGORITHMS_H
#define MATH_ALGORITHMS_H

#include <stdlib.h>
#include <math.h>

#ifndef M_SQRT2
	#define M_SQRT2     1.41421356237309504880
#endif
#ifndef M_PI
	#define M_PI		3.14159265358979323846
#endif

double antierf(double x, double eps = 1.e-8){
    if (fabs(x) > 1) return NAN;
    if (x == 1) return INFINITY;
    if (x == -1) return -INFINITY;
    double y = 0.;
    double dy = 0.;
    // Число циклов < 32 при x = 0.999999999999 y ~ 5
    do {
        y += dy;
        dy = (x - erf(y))/(2./sqrt(M_PI)*exp(-y*y));
    } while (fabs(dy) > eps*fabs(y));
    return y;
}

double drand(double l = 1.){
    return l*rand()/RAND_MAX;
}

double nrand(double m = 0., double sigma = 1.){
    return m + M_SQRT2*sigma*antierf(2*drand() - 1);
}

#endif // MATH_ALGORITHMS_H
