#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <time.h>

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
          )
{
    double r, dz;
    switch (Method){
    case 1:
        for (int i = 1; i<N-1; i++){
            dz = z[i] - z[i-1];
            r = dt*dt/dz/dz*cdz2(z[i],t);
            if ((r<0)||(r>1)) {
                printf("Плохое число Куранта r = %f", r);
                return;
            }
            u[i] = r*(up[i+1] + up[i-1]) + 2*(1 - r)*up[i] + source(z[i],t)*dt*dt - upp[i];
        }
        dz = z[1] - z[0];
        u[0] = (b1 - dc1/dz*u[1])/(c1 - dc1/dz);
        dz = z[N - 1] - z[N - 2];
        u[N - 1] = (b2 + dc2/dz*u[N - 2])/(c2 + dc2/dz);
        break;
    case 2:
        double A[N], B[N], C[N];
        dz = z[1] - z[0];
        A[0] = 0;
        B[0] = c1 - dc1/dz;
        C[0] = dc1/dz;
        u[0] = b1;
        for (int i = 1; i<N-1; i++){
            dz = z[i] - z[i-1];
            r = dt*dt/dz/dz*cdz2(z[i],t);
            A[i] = -r;
            B[i] = 1 + 2*r;
            C[i] = -r;
            u[i] = 2*up[i] - upp[i] + dt*dt*source(z[i], t);
        }
        dz = z[N - 1] - z[N - 2];
        A[N-1] = -dc2/dz;
        B[N-1] = c2 + dc2/dz;
        C[N-1] = 0;
        u[N-1] = b2;
        for (int i = 0; i<N-1; i++){
            u[i+1] -= u[i]*A[i+1]/B[i];
            B[i+1] -= C[i]*A[i+1]/B[i];
        }
        u[N-1] *= 1./B[N-1];
        for (int i = N - 1; i>0; i--){
            u[i-1] = (u[i-1] - u[i]*C[i-1])/B[i-1];
        }
        break;
    }
};

double ident(double z, double t){
    return 1.;
};

double null(double z, double t){
    return 0.;
};

void he(double *u, int N, int M){
    double lmax = 1., tmax = 1.;
    double z[N];
    double dz = lmax/(N - 1), dt = tmax/(M - 2), t = 0.;
    for (int i = 0; i<N; i++){
        z[i] = i*dz;
        u[1*N + i] = 1 - 4*(z[i] - 0.5)*(z[i] - 0.5);
        //u[1*N + i] = sin(M_PI*z[i]);
        //u[1*N + i] = (fabs(z[i] - 0.5)<0.5) ? 1 : 0;
        u[0*N + i] = u[1*N + i] - 0*dt;
    }
    for (int j = 2; j<M; j++){
        t += dt;
        pdehe(2, &u[j*N], &u[(j-1)*N], &u[(j-2)*N], N, z, t, dt, ident, null, 1., 0., 1., 0., 0., 0.);
    }
}

typedef std::complex<double> clex;

void ft(clex *x, clex *y, int N){
    for (int i = 0; i<N; i++){
        y[i] = 0;
        for (int j = 0; j<N; j++){
            y[i] += x[j]*exp(clex(0, 2*M_PI*j*i/N));
        }
        //printf("%f %f\n", real(y[i]), imag(y[i]));
    }
};

void ft_inverse(clex *y, clex *x, int N){
    for (int i = 0; i<N; i++){
        x[i] = 0;
        for (int j = 0; j<N; j++){
            x[i] += 1./N*y[j]*exp(clex(0, -2*M_PI*j*i/N));
        }
        //printf("%f %f\n", real(y[i]), imag(y[i]));
    }
};

clex gauss(clex x, double av, double sigma){
    return 1/sqrt(2*M_PI)*exp(-pow(x - av, 2)/2./sigma/sigma);
}

clex impulse(clex x, double av, double length){
    return (fabs(real(x) - av)<length/2) ? 1 : 0;
}

double random(){
    return 1.0*rand()/RAND_MAX;
}

void fourie(){
    int N = 100, M = N;
    clex x[N], y[N], xi[N];
    double *u;
    u = new double[M*N];
    srand (time(NULL));
    double z[N], omega[N], dt = 1./N;
    int number = 11;
    he(u, N, M);
    while (number<13){
        for (int i = 0; i<N; i++){
            z[i] = i*dt;
            omega[i] = 1.*i/dt/N;
            switch (number){
            case 1:
                x[i] = cos(4*z[i]*M_PI*2);
                break;
            case 2:
                x[i] = gauss(clex(z[i]), 0.5, 0.1);
                break;
            case 3:
                x[i] = gauss(clex(z[i]), 0.5, 0.03);
                break;
            case 4:
                x[i] = impulse(clex(z[i]), 0.5, 0.03);
                break;
            case 5:
                x[i] = impulse(clex(z[i]), 0.5, 0.3);
                break;
            case 6:
                x[i] = 0.1*random() + impulse(clex(z[i]), 0.5, 0.3);
                break;
            case 7:
                x[i] = 0.1*random() + gauss(clex(z[i]), 0.5, 0.1);
                break;
            case 8:
                x[i] = 0.1*random() + cos(4*z[i]*M_PI*2);
                break;
            case 9:
                x[i] = cos(4*z[i]*M_PI*2) + cos(8*z[i]*M_PI*2);
                break;
            case 10:
                x[i] = (z[i] < 0.5) ? cos(4*z[i]*M_PI*2) : cos(8*z[i]*M_PI*2);
                break;
            case 11:
                x[i] = u[(M-1)*N + i];
                break;
            case 12:
                x[i] = u[i*N + N/2];
                break;
            case 13:
                x[i] = u[(M-1)*N + i];
                break;
            case 14:
                x[i] = u[i*N + N/2];
                break;
            case 15:
                x[i] = u[(M-1)*N + i];
                break;
            case 16:
                x[i] = u[i*N + N/2];
                break;

            }
            //gauss(clex(z[i]), 0.5, 0.1);//sin(z[i]*M_PI*2);//
        }
        ft(x, y, N);
        ft_inverse(y, xi, N);

        FILE *fd;
        fd = fopen("lab4.dat", "w");
        for (int i = 0; i<N; i++){
            fprintf(fd, "%f %f %f %f %f %f %f %f\n", z[i], real(x[i]), omega[i], real(y[i]), imag(y[i]), abs(y[i]), arg(y[i]), real(xi[i]));
        }
        fclose(fd);

        fd = fopen("lab4.gp", "w");
        fprintf(fd, "set term png size 800,800 enhanced font 'Verdana,10'\n");
        fprintf(fd, "set output 'f/lab4-%d-f.png'\n", number);
        fprintf(fd, "set multiplot layout 2,2 rowsfirst\n");
        fprintf(fd, "unset key \n");
        fprintf(fd, "set grid \n");
        fprintf(fd, "plot 'lab4.dat' using 1:2 w l lw 2, '' using 1:8 w l lw 2\n");
        fprintf(fd, "plot 'lab4.dat' using 3:4 w impulses lw 2, '' using 3:5 w impulses lw 2 \n");
        fprintf(fd, "plot 'lab4.dat' using 3:6 w impulses lw 2\n");
        fprintf(fd, "plot 'lab4.dat' using 3:7 w impulses lw 2\n");
        fclose(fd);

        fd = popen("lab4.gp", "r");
        pclose(fd);
        number++;
    }
}

clex morle(clex x, double alpha, double k0){
    return exp(-x*x/alpha/alpha)*(exp(clex(0., k0)*x));
}

void wavelet_morle(double *t, clex *x, int N,  clex *y, double *a, int Na, double *b, int Nb, double alpha, double k0){
    // y размера Na*Nb
    clex n, w;
    for (int i = 0; i<Na; i++){
        for (int j = 0; j<Nb; j++){
            w = 0;
            n = 0;
            for (int k = 0; k<N; k++){
                w += x[k]*morle((t[k] - b[j])/a[i], alpha, k0);
                n += exp(-1./alpha/alpha*pow((t[k] - b[j])/a[i], 2));
            }
            y[i*Nb + j] = w/n;
        }
    }
};

void wavelet(){
    int N = 100;
    int Na = 100;
    int Nb = 200;
    int M = N;
    double dt = 1./(N - 1);

    double alpha = sqrt(2);
    double k0 = M_PI*sqrt(2/log(2));

    double amin = 2*dt/alpha;
    double amax = 1./alpha;
    double bmin = 0;
    double bmax = (N - 1)*dt;
    double a[Na], b[Nb], z[N];
    double *u;
    u = new double[M*N];
    for (int i = 0; i<Na; i++) a[i] = amin + i*(amax - amin)/(Na - 1);
    for (int i = 0; i<Nb; i++) b[i] = bmin + i*(bmax - bmin)/(Nb - 1);
    printf("%e \n", a[0]);
    clex x[N], *y;
    y = new clex[Na*Nb];
    he(u, N, M);

    int number = 18;
    while (number<20){
        for (int i = 0; i<N; i++){
            z[i] = i*dt;
            switch (number){
            case 1:
                x[i] = cos(4*z[i]*M_PI*2);
                break;
            case 2:
                x[i] = gauss(clex(z[i]), 0.5, 0.1);
                break;
            case 3:
                x[i] = gauss(clex(z[i]), 0.5, 0.03);
                break;
            case 4:
                x[i] = impulse(clex(z[i]), 0.5, 0.03);
                break;
            case 5:
                x[i] = impulse(clex(z[i]), 0.5, 0.3);
                break;
            case 6:
                x[i] = 0.1*random() + impulse(clex(z[i]), 0.5, 0.3);
                break;
            case 7:
                x[i] = 0.1*random() + gauss(clex(z[i]), 0.5, 0.1);
                break;
            case 8:
                x[i] = 0.1*random() + cos(4*z[i]*M_PI*2);
                break;
            case 9:
                x[i] = cos(2*z[i]*M_PI*2) + cos(4*z[i]*M_PI*2);
                break;
            case 10:
                x[i] = (z[i] < 0.5) ? cos(2*z[i]*M_PI*2) : cos(4*z[i]*M_PI*2);
                break;
            case 11:
                x[i] = (z[i] < 0.5) ? sin((z[i] - 0.5)*M_PI*2) : sin(3*(z[i] - 0.5)*M_PI*2);
                break;
            case 12:
                x[i] = morle(z[i], alpha, k0) + morle(z[i], alpha, 2*k0);
                break;
            case 13:
                x[i] = morle(z[i], alpha, k0) + morle(z[i] - 0.5, alpha, k0);
                break;
            case 14:
                x[i] = u[(M-1)*N + i];
                break;
            case 15:
                x[i] = u[i*N + N/2];
                break;
            case 16:
                x[i] = u[(M-1)*N + i];
                break;
            case 17:
                x[i] = u[i*N + N/2];
                break;
            case 18:
                x[i] = u[(M-1)*N + i];
                break;
            case 19:
                x[i] = u[i*N + N/2];
                break;
            }
        }
        wavelet_morle(z, x, N,  y, a, Na, b, Nb, alpha, k0);

        FILE *fd;
        fd = fopen("lab4.dat", "w");
        for (int i = 0; i<Na; i++){
            for (int j = 0; j<Nb; j++){
                fprintf(fd, "%f %f %f\n", b[j], k0/a[i], abs(y[i*Nb + j]));
            }
            fprintf(fd, "\n");
        }
        fclose(fd);

        fd = fopen("lab4.gp", "w");
        fprintf(fd, "set term png size 800,800 enhanced font 'Verdana,10'\n");
        fprintf(fd, "set output 'w2/lab4-%d-w.png'\n", number);
        //fprintf(fd, "f(x,y)=sin(1.3*x)*cos(.9*y)+cos(.8*x)*sin(1.9*y)+cos(y*.2*x)\n");
        fprintf(fd, "set xrange [0:1]\n");
        //fprintf(fd, "set yrange [0:%e]\n", k0);
        fprintf(fd, "set isosample 250, 250\n");
        fprintf(fd, "set table 'test.dat'\n");
        fprintf(fd, "splot 'lab4.dat' using 1:2:3\n");
        fprintf(fd, "unset table\n");
        fprintf(fd, "set contour base\n");
        fprintf(fd, "set cntrparam\n");
        fprintf(fd, "unset surface\n");
        fprintf(fd, "set table 'cont.dat'\n");
        fprintf(fd, "splot 'test.dat'\n");
        fprintf(fd, "unset table\n");
        fprintf(fd, "reset\n");
        //fprintf(fd, "set xrange [-5:5]\n");
        //fprintf(fd, "set yrange [-5:5]\n");
        fprintf(fd, "unset key\n");
        fprintf(fd, "unset clabel\n");
        fprintf(fd, "set palette rgbformulae 33,13,10\n");
        fprintf(fd, "p 'test.dat' with image, 'cont.dat' w l lt -1 lw 1.5\n");
        fclose(fd);

        fd = popen("lab4.gp", "r");
        pclose(fd);
        number++;
    }
}

int main(){
    fourie();
    //wavelet();
    return 0;
}
