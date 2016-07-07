#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
	#define M_PI		3.14159265358979323846
#endif

///////////////////////////////////////////////////////////////////
// Math functions
///////////////////////////////////////////////////////////////////

double arsh(double x) {
	return log(x+sqrt(x*x+1));
}

double arth(double x) {
	return 1.0/2*log((1+x)/(1-x));
}

double random(const double &xmin, const double &xmax){
    return xmin + (xmax - xmin)*1.*rand()/RAND_MAX;
}

template <typename T>
T min(T a, T b) {
    return (a > b) ? b : a;
}

template <typename T>
T max(T a, T b) {
    return (a > b) ? a : b;
}
