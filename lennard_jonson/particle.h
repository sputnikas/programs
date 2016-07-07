#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdio.h>
#include <math.h>

const double CL = 2.99792458E8;
const double CL2 = 8.9875517873681764E16;

struct Particle{
	Vector r;
	Vector v;
	double m;
	double R;
};

double energy(const Particle &p){
	return p.m*CL2/sqrt(1 - dot(p.v, p.v)/CL2);
}

double energy_p(const Vector &p, const double &m){
	return sqrt(dot(p, p)*CL2 + pow(m*CL2, 2));
}

double energy(const Vector &v, const double &m){
	return m*CL2/sqrt(1 - dot(v, v)/CL2);
}

Vector impulse(const Particle &p){
	return p.m/sqrt(1 - dot(p.v, p.v)/CL2)*p.v;
}

Vector impulse(const Vector &v, const double &m){
	return m/sqrt(1 - dot(v, v)/CL2)*v;
}

Vector velocity(const Vector &p, const double &m){
    return p*CL2/energy_p(p, m);
}

Vector a(const Vector &F, const Vector &v, const double &m){
    if (!(sqrt(1 - dot(v, v)/CL2) == sqrt(1 - dot(v, v)/CL2))){
        return {0., 0., 0.};
    }
    return (F - v*dot(v, F)/CL2)*sqrt(1 - dot(v, v)/CL2)/m;
}

void to_console(Particle p){
    printf("r = {%2.4e, %2.4e, %2.4e} ", p.r.x, p.r.y, p.r.z);
    printf("v = {%2.4e, %2.4e, %2.4e} ", p.v.x, p.v.y, p.v.z);
    printf("m = %2.4e ", p.m);
    printf("\n");
}


#endif // PARTICLE_H
