#ifndef VECTOR_H
#define VECTOR_H

#include <stdio.h>
#include <math.h>

struct Vector{
	double x;
	double y;
	double z;
	Vector(): x(0), y(0), z(0) {};
	Vector(double x, double y, double z) : x(x), y(y), z(z) {};
};

Vector operator * (const Vector &a, const Vector &b){
	return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

Vector operator + (const Vector &a, const Vector &b){
	return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vector& operator += (Vector &a, const Vector &b){
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}

Vector& operator -= (Vector &a, const Vector &b){
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}

Vector operator - (const Vector &a, const Vector &b){
	return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vector operator - (const Vector &a){
	return {- a.x, - a.y, - a.z};
}

Vector operator * (const double &p, const Vector &a){
	return {p*a.x, p*a.y, p*a.z};
}

Vector operator * (const Vector &a, const double &p){
	return {p*a.x, p*a.y, p*a.z};
}

Vector operator / (const Vector &a, const double &p){
	return {a.x/p, a.y/p, a.z/p};
}

double dot(const Vector &a, const Vector &b){
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

double norm(const Vector &a){
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

void to_console(Vector v){
    printf("{%2.4e, %2.4e, %2.4e}", v.x, v.y, v.z);
}


#endif // VECTOR_H
