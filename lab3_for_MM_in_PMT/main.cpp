#include "main.h"

#define PI 3.14169
double epsilon = 1.2e-2;

double null(double z){
    return 0;
}

double delta(double z){
    return ((z>-epsilon/2)&&(z<epsilon/2)) ? 1/epsilon : 0;
}

double s1(double z, double t){
    return delta(z - 0.5);
}

double cd1(double z, double t){
    return 0.1*sin((z - 0.5)*M_PI*6) + 1;
}

double ident(double z, double t = 0){
    return 1;
}

double g1(double z, double t){
    return (sin(PI*t)+1)*0.5;
}

double g2(double z, double t){
    return 1 + (sin(0.5*PI*t)+1)*0.05;
}

double f1(double z){
    return sin(PI*z);
}

double f2(double z){
    return 1 - sin(PI*z);
}

double f3(double z){
    return 4*(-(z - 0.5)*(z - 0.5) + 0.25);
}

double f5(double z){
    return 0.5*(1 + sin(PI*z));
}

double f6(double z){
    return sin(3*PI*z);
}

double f7(double z){
    return (z<0.5) ? z : 1 - z;
}

double f8(double z, double z0, double a){
    return (fabs(z - z0)<a/2) ? 1 : 0;
}

double f9(double x){
    return 1;
}

double f10(double x){
    return 2;
}

double es1(double x, double y, double u){
    return delta(x - 0.5)*delta(y - 0.5);
}

double es2(double x, double y, double u){
    return delta(x - 0.3)*delta(y - 0.3) - delta(x - 0.7)*delta(y - 0.7);
}

int main(int argv, char **argc){
    Diffuse();
    return 0;
};

void Elleptic(){
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    //printf("%d \n", argv);

    int N = 100;
    int M = 100;
    double *u, x[N], y[M];
    u = new double[N*M];
    double x0 = 0, xmax = 1.;
    for (int i = 0; i<N; i++) x[i] = x0 + i*(xmax - x0)/(N - 1);
    double y0 = 0, ymax = 1.;
    for (int i = 0; i<M; i++) y[i] = y0 + i*(ymax - y0)/(M - 1);

    FILE *gp;
    gp = fopen("e5.gpl", "w");
    //fprintf(gp, "reset\n");
    fprintf(gp, "set dummy u,v\n");
    fprintf(gp, "set samples 51, 51\n");
    fprintf(gp, "set isosamples 21, 21\n");
    //fprintf(gp, "set palette rgbformulae 20, 25, 25\n");
    fprintf(gp, "set contour\n");
    fprintf(gp, "set xrange [ 0.00000 : 1.00000 ] noreverse nowriteback\n");
    fprintf(gp, "set yrange [ 0.00000 : 1.00000 ] noreverse nowriteback\n");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set pm3d at s\n");
    fprintf(gp, "set hidden3d\n");
    fprintf(gp, "splot '-' using 1:2:3 with lines \n");
    int i1 = N/3, j1 = M/3;
    int i2 = 2*N/3, j2 = 2*M/3;
    double u1 = 2, u2 = 4;
    pdeee(u, N, M, x, y, 0.9, 0.0000003,
          i1, j1, u1,
          i2, j2, u2,
          0, 1, null,
          0, 1, null,
          0, 1, null,
          0, 1, null
          );
    double current_h = 0, dx, dy;
    if (i2 - i1 > 3){
        int iav = (i1 + i2)/2;
        dx = x[iav + 1] - x[iav];
        for (int j = 1; j<M; j++){
            dy = y[j] - y[j-1];
            current_h += (u[(iav + 1)*M + j] - u[iav*M + j])/dx*dy;
        }
    }
    fprintf(gp, "# I/h = %f U = %f R*h = %f\n", current_h, u2 - u1, (u2 - u1)/current_h);
	for (int j = 0; j<M-1; j++){
        for (int i = 0; i<N; i++) fprintf(gp, "%f %f %f\n", x[i], y[j], u[i*M + j]);
        fprintf(gp, "\n");
	}
	//fclose(gp);
    //gp = popen("lab3.gp", "r");
    //pclose(gp);
}

void Diffuse(){
    int N = 1000, M = 1000;
    int k = 10;
    double *u, *up, *upp, z[N], t[M], *temp;
    u = new double[N];
    up = new double[N];
    double z0 = 0, zmax = 1.;
    for (int i = 0; i<N; i++) z[i] = z0 + i*(zmax - z0)/(N - 1);
    double t0 = 0, tmax = 10.;
    for (int i = 0; i<M; i++) t[i] = t0 + i*(tmax - t0)/(M - 1);

    for (int i = 0; i<N; i++) {
        u[i] = null(z[i]);
        //up[i] = u[i] - null(z[i])*(t[1] - t[0]);
    }
    FILE *gp, *gnuplot,*tmp;
    gp = fopen("lab3.gp", "w");
    fprintf(gp, "set term gif ");
    fprintf(gp, "animate ");
    //fprintf(gp, "optimize ");
    fprintf(gp, "delay 1 ");
    fprintf(gp, "size 600, 600 ");
    fprintf(gp, "background \"#ffffff\" ");
    //fprintf(gp, "crop ");
    fprintf(gp, "\n");
    fprintf(gp, "set output 'PE11.gif' \n ");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set xrange[0:1]\n");
    fprintf(gp, "set yrange[0:4.2]\n");
    fprintf(gp, "set grid x y\n");
	for (int j = 0; j<M-1; j++){
        if (j%k == 0) {
            if (j == 0)
                fprintf(gp, "plot '-' using 1:2 with lines lw 2 \n");
            else
                fprintf(gp, "plot '-' using 1:2 with lines lw 2  \n");
            for (int i = 0; i<N; i++) fprintf(gp, "%f %f\n", z[i], u[i]);
            fprintf(gp, "end\n");
            //fprintf(gp, "replot\n");
        }
        temp = up;
        up = u;
        u = temp;
        //for (int i = 0; i<N; i++) printf("%f ", uprev[i]);
        //system("pause");
        pdepe(2, u, up, N, z, t[j], t[j+1] - t[j], cd1 , ident, 1, 0, 1, 0, 1, 1);
        //pdehe(2, u, up, upp, N, z, t[j], t[j+1] - t[j],  ident, null, 1, 0, 1, 0, 0, 0);
	}
	fclose(gp);
    gp = popen("lab3.gp", "r");
    pclose(gp);

}
