#ifndef MAIN_H
#define MAIN_H

#include <stdio.h>
#include <math.h>

#include "vector.h"
#include "particle.h"
#include "math_algorithms.h"

#ifndef M_PI
	#define M_PI		3.14159265358979323846
#endif

#define N_PARAMS 10
#define N_COUNTS 500

#define ERR_OFILE "Unable to open file %s\n"
#define ERR_CFILE "Unable to create file %s\n"
#define CREATE_FILE "File %s has been created\n"
#define ERR_FILEEXIST "Check, that file %s exists (in the directory from which application was run)\n"
#define NAME_FILE "Name of file <%s>\n"
#define SECS "Count time: %.1f secs\n"

const double K = 1.38E-23;

int collide_elastic(Particle &p1, Particle &p2, double dt){
    Vector r12 = p1.r - p2.r;
    Vector v12 = p1.v - p2.v;
    double nr12 = norm(r12);
    double nv12 = norm(v12);
    double R = p1.R + p2.R;
    double cos_theta = dot(v12, r12)/nr12/nv12;
    double sin_theta = sqrt(1 - cos_theta*cos_theta);
    double tau = nr12/nv12*(- cos_theta - sqrt(pow(R/nr12, 2) - pow(sin_theta, 2)));

    if ((1 >= R/nr12)&&(cos_theta <= 0)&&(R/nr12 >= fabs(sin_theta))&&(tau < dt)){
        Vector n = (r12 + v12*tau)/norm(r12 + v12*tau);
        Vector i1 = impulse(p1);
        Vector i2 = impulse(p2);
        double i1n = dot(i1, n);
        double i2n = dot(i2, n);
        double in = i1n + i2n;
        double e1 = energy(p1)/CL;
        double e2 = energy(p2)/CL;
        double e = e1 + e2;
        Vector i12 = 2*(i1n*e2 - i2n*e1)*e/(e*e - in*in)*n;
        //p1.r = p1.r + p1.v*tau;
        //p2.r = p2.r + p2.v*tau;
        p1.v = velocity(i1 - i12, p1.m);
        p2.v = velocity(i2 + i12, p2.m);
        //p1.r = p1.r + p1.v*(dt - tau);
        //p2.r = p2.r + p2.v*(dt - tau);
        return 1;
    }
    return 0;
}

///////////////////////////////////////////////////////////////////
// Struct of file
///////////////////////////////////////////////////////////////////

struct Params{
	double dt;              // Шаг по времени
	int method_num;         // 1 - метод Рунге-Кутты 4-го порядка; 1 - метод Рунге-Кутты 4-го порядка по импульсам
	int n_particle;         // Число частиц
	int in_field;           // Учитывать или не учитывать взаимодействие? 0 - не учитывать; 1 - Леннард-Джонсон; 2 - Ньютонова гравитация,
	int ex_field;

	double T_e;
	double T_i;
	double m_e;
	double m;
	double r;

	int type_space;         // 0 - куб, 1 - сфера.
	double a;               // Сторона куба, радиус сферы

    double in[N_PARAMS];
    double ex[N_PARAMS];
};

void to_console(Params h){
    printf("h.dt         = %e\n", h.dt         );
    printf("h.method_num = %d\n", h.method_num );
    printf("h.n_particle = %d\n", h.n_particle );
    printf("h.in_field   = %d\n", h.in_field   );
    printf("h.ex_field   = %d\n", h.ex_field   );
    printf("h.T_e        = %e\n", h.T_e        );
    printf("h.T_i        = %e\n", h.T_i        );
    printf("h.m_e        = %e\n", h.m_e        );
    printf("h.m          = %e\n", h.m          );
    printf("h.r          = %e\n", h.r          );
    printf("h.type_space = %d\n", h.type_space );
    printf("h.a          = %e\n", h.a          );
}

extern Params h;
extern Particle *particles[N_COUNTS];

int load_ini(const char *filename){
    struct Ini{
        char key[150];
        char val[150];
    };
    int n = 0;
    Ini in[100];

    FILE *f;
	if ((f = fopen(filename, "r")) == NULL) {
		printf(ERR_OFILE, filename);
		return 0;
	}

    char buf[300];
	char parts[2][150];
	char *buf2;

	while (fgets(buf, sizeof(buf), f) != NULL){
        if (strspn((const char*)buf, "; ") == strspn((const char*)buf, " ")) {
            if (sscanf(buf, "%s = %[^;\n]", parts[0], parts[1]) == 2) {
                strcpy(in[n].key, parts[0]);
                strcpy(in[n].val, parts[1]);
                n++;
            }
        }
	}

    for (int i = 0; i<n; i++){
        if (strcmp("h.dt", in[i].key) == 0) sscanf(in[i].val, "%le", &h.dt);
        if (strcmp("h.method_num", in[i].key) == 0) sscanf(in[i].val, "%d", &h.method_num);
        if (strcmp("h.n_particle", in[i].key) == 0) sscanf(in[i].val, "%d", &h.n_particle);
        if (strcmp("h.in_field", in[i].key) == 0) sscanf(in[i].val, "%d", &h.in_field);
        if (strcmp("h.ex_field", in[i].key) == 0) sscanf(in[i].val, "%d", &h.ex_field);
        if (strcmp("h.type_space", in[i].key) == 0) sscanf(in[i].val, "%d", &h.type_space);
        if (strcmp("h.a", in[i].key) == 0) sscanf(in[i].val, "%le", &h.a);
        if (strcmp("h.T_e", in[i].key) == 0) sscanf(in[i].val, "%le", &h.T_e);
        if (strcmp("h.T_i", in[i].key) == 0) sscanf(in[i].val, "%le", &h.T_i);
        if (strcmp("h.m", in[i].key) == 0) sscanf(in[i].val, "%le", &h.m);
        if (strcmp("h.r", in[i].key) == 0) sscanf(in[i].val, "%le", &h.r);
        if (strcmp("h.m_e", in[i].key) == 0) sscanf(in[i].val, "%le", &h.m_e);
        if (strcmp("in", in[i].key) == 0) {
            int k = 0;
            buf2 = strtok(in[i].val, ", \t");
            do{
                if (k<N_PARAMS) sscanf(buf2, "%le", &h.in[k]);
                k++;
                buf2 = strtok(NULL, ", \t");
            } while ((buf2)&&(k<N_PARAMS));
        }
        if (strcmp("ex", in[i].key) == 0) {
            int k = 0;
            buf2 = strtok(in[i].val, ", \t");
            do{
                if (k<N_PARAMS) sscanf(buf2, "%le", &h.ex[k]);
                k++;
                buf2 = strtok(NULL, ", \t");
            } while ((buf2)&&(k<N_PARAMS));
        }
    }

    fclose(f);
    return 0;
}

void runge( Particle &next,
            const Particle &p,
            const double &t,
            const double &dt,
            const Vector &F,
            Vector (*force)(Vector, Vector, double, double))
{
	Vector k1v,  k2v, k3v, k4v;
	Vector k1r,  k2r, k3r, k4r;
	k1v = dt*a(force(p.r, p.v, t, p.m) + F, p.v, p.m);
	k1r = dt*p.v;
	k2v = dt*a(force(p.r + k1r/2, p.v + k1v/2, t + dt/2, p.m) + F, p.v + k1v/2, p.m);
	k2r = dt*(p.v + k1v/2);
	k3v = dt*a(force(p.r + k2r/2, p.v + k2v/2, t + dt/2, p.m) + F, p.v + k2v/2, p.m);
	k3r = dt*(p.v + k2v/2);
	k4v = dt*a(force(p.r + k3r, p.v + k3v, t + dt, p.m) + F, p.v + k3v, p.m);
	k4r = dt*(p.v + k3v);
	to_console(k1r);
	next.v = p.v + 1.0/6*(k1v+2*k2v+2*k3v+k4v);
	next.r = p.r + 1.0/6*(k1r+2*k2r+2*k3r+k4r);
};

void runge_impulse( Particle &next,
                    const Particle &p,
                    const double &t,
                    const double &dt,
                    const Vector &F,
                    Vector (*force)(Vector, Vector, double, double))
{
	Vector k1p, k2p, k3p, k4p;
	Vector k1r,  k2r, k3r, k4r;
	Vector imp = impulse(p);
	k1p = dt*(force(p.r, velocity(imp, p.m), t, p.m) + F);
	k1r = dt*velocity(imp, p.m);
	k2p = dt*(force(p.r + k1r/2, velocity(imp + k1p/2, p.m), t + dt/2, p.m) + F);
	k2r = dt*velocity(imp + k1p/2, p.m);
	k3p = dt*(force(p.r + k2r/2, velocity(imp + k2p/2, p.m), t + dt/2, p.m) + F);
	k3r = dt*velocity(imp + k2p/2, p.m);
	k4p = dt*(force(p.r + k3r, velocity(imp + k3p, p.m), t + dt, p.m) + F);
	k4r = dt*velocity(imp + k3p, p.m);
	next.v = velocity(impulse(p.v, p.m) + 1.0/6*(k1p+2*k2p+2*k3p+k4p), p.m);
	next.r = p.r + 1.0/6*(k1r+2*k2r+2*k3r+k4r);
};

Vector force(Vector r, Vector v, double t, double m){
    switch (h.ex_field){
    default:
        return {0., 0., 0.};
        break;
    case 1:
        return m*Vector(0., 0., -h.ex[0]);
        break;
    }
}

Vector force(const Particle &p1, const Particle &p2){
    Vector r = p1.r - p2.r;
    double R = norm(r);
    switch (h.in_field){
    default:
        return {0., 0., 0.};
        break;
    case 1:
        return r*(- h.in[0]*h.in[2]/pow(R, h.in[0] + 2) + h.in[1]*h.in[3]/pow(R, h.in[1] + 2));
        break;
    case 2:
        return h.in[0]*r*p1.m*p2.m/pow(R, h.in[1]);
        break;
    }
}

void collide_wall(Particle &p1, double dt, Vector n){
    double vq = sqrt(K*h.T_e/2/h.n_particle/h.m_e);
    double pq = h.m_e*vq*sqrt(2 + vq*vq/CL2);

    Vector i1 = impulse(p1);
    Vector i2 = Vector(nrand(0., pq), nrand(0., pq), nrand(0., pq));
    double i1n = dot(i1, n);
    double i2n = fabs(nrand(0., pq));
    while (i2n < fabs(i1n))
        i2n = fabs(nrand(0., pq));
    i2 = i2 - dot(i2, n)*n + i2n*n;
    double in = i1n + i2n;
    double e1 = energy(p1)/CL;
    double e2 = energy_p(i2, h.m_e)/CL;
    double e = e1 + e2;
    Vector i12 = 2*(i1n*e2 - i2n*e1)*e/(e*e - in*in)*n;
    //p1.r = p1.r + p1.v*tau;
    p1.v = velocity(i1 - i12, p1.m);
    //p1.r = p1.r + p1.v*(dt - tau);
}

int calculate(int nt){
    double t = h.dt*nt;
    double dt = h.dt;
    Vector F, r;
    if (nt < N_COUNTS - 1){
        for (int i = 0; i<h.n_particle; i++){
            particles[nt + 1][i] = particles[nt][i];
            switch (h.type_space){
            case 0:
                r = particles[nt + 1][i].r + particles[nt + 1][i].v*h.dt;
                if (r.x > h.a/2)
                    collide_wall(particles[nt + 1][i], h.dt, Vector(-1., 0., 0.));
                if (r.x < - h.a/2)
                    collide_wall(particles[nt + 1][i], h.dt, Vector(1., 0., 0.));
                if (r.y > h.a/2)
                    collide_wall(particles[nt + 1][i], h.dt, Vector(0., -1., 0.));
                if (r.y < - h.a/2)
                    collide_wall(particles[nt + 1][i], h.dt, Vector(0., 1., 0.));
                if (r.z > h.a/2)
                    collide_wall(particles[nt + 1][i], h.dt, Vector(0., 0., -1.));
                if (r.z < - h.a/2)
                    collide_wall(particles[nt + 1][i], h.dt, Vector(0., 0., 1.));
                break;
            case 1:
                r = particles[nt + 1][i].r + particles[nt + 1][i].v*h.dt;
                if (norm(r) > h.a)
                    collide_wall(particles[nt + 1][i], h.dt, r/norm(r));
                break;
            default:
                break;
            }
            F = {0., 0., 0.};
            for (int j = 0; (j<h.n_particle)&&(j!=i); j++){
                if (collide_elastic(particles[nt + 1][i], particles[nt + 1][j], h.dt) == 0)
                    F = F + force(particles[nt + 1][i], particles[nt + 1][j]);
            }
            switch (h.method_num){
            default:
                runge(particles[nt + 1][i], particles[nt + 1][i], t, dt, F, force);
                break;
            case 2:
                runge_impulse(particles[nt + 1][i], particles[nt + 1][i], t, dt, F, force);
                break;
            }
        }
    }
    return 1;
}

void to_file(const char *name){
    FILE *f;
    f = fopen(name, "w");

    for (int i = 0; i<N_COUNTS; i++){
        for (int j = 0; j<h.n_particle; j++){
            fprintf(f, "(%e, %e, %e, %e, %e, %e)\v",
                    particles[i][j].r.x,
                    particles[i][j].r.y,
                    particles[i][j].r.z,
                    particles[i][j].v.x,
                    particles[i][j].v.y,
                    particles[i][j].v.z);
        }
        fprintf(f, "\n");
    }

    fclose(f);
}

void to_bfile(const char *name){
    FILE *f;
    f = fopen(name, "wb");

    fwrite(&h, sizeof(h), 1, f);
    for (int i = 0; i<N_COUNTS; i++){
        fwrite(&particles[i][0], sizeof(Particle), h.n_particle, f);
    }

    fclose(f);
}

void from_bfile(const char *name){
    FILE *f;
    f = fopen(name, "rb");

    fread(&h, sizeof(h), 1, f);
    for (int i = 0; i<N_COUNTS; i++){
        fread(&particles[i][0], sizeof(Particle), h.n_particle, f);
    }

    fclose(f);
}

int init(const char* name){
    load_ini(name);
    to_console(h);
    for (int i = 0; i<N_COUNTS; i++) particles[i] = new Particle[h.n_particle];
    int a = 0, b = 0, c = 0;
    int cub = (int) trunc(pow(h.n_particle, 1./3)) + 1;
    double cub_a = cub*2*h.r;
    double cub_da = cub_a / (cub - 1);
    double vq = sqrt(K*h.T_i/2/h.n_particle/h.m);
    double pq = h.m*vq*sqrt(2 + vq*vq/CL2);
    for (int i = 0; i<h.n_particle; i++){
        particles[0][i].R = h.r;
        particles[0][i].m = h.m;
        particles[0][i].r = {-cub_a/2 + cub_da*a, -cub_a/2 + cub_da*b, -cub_a/2 + cub_da*c};
        a++;
        if (a >= cub){
            a = 0;
            b++;
            if (b >= cub){
                b = 0;
                c++;
            }
        }
        particles[0][i].v = velocity(Vector(nrand(0., pq), nrand(0., pq), nrand(0., pq)), h.m);
    }
    to_file("test1.particles");
    return 0;
}


int free(){
    for(int i = 0; i<N_COUNTS; i++) delete [] particles[i];
    return 0;
}

#endif // MAIN_H
