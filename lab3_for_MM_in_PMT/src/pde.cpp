#include "pde.h"

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
          )
{
    double r, dz;
    switch (Method){
    case 1:
        for (int i = 1; i<N-1; i++){
            dz = z[i] - z[i-1];
            r = dt/dz/dz*cdz2(z[i],t);
            if ((r<0)||(r>1)) {
                printf("Плохое число Куранта r = %f", r);
                return;
            }
            u[i] = r*(up[i+1] + up[i-1]) + (1 - 2*r)*up[i] + source(z[i],t)*dt;
            //printf("%f ", up[i]);
        }
        //printf("\n");
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
            dz = z[i+1] - z[i];
            r = dt/dz/dz*cdz2(z[i],t);
            A[i] = r;
            B[i] = - 1 - 2*r;
            C[i] = r;
            u[i] = - up[i] - dt*source(z[i], t);
        }
        //printf("%f ", r);
        A[N-1] = - dc2/dz;
        B[N-1] = c2 + dc2/dz;
        C[N-1] = 0;
        u[N-1] = b2;
        //printf("%f %f ", up[N-2], u[N-2]);
        for (int i = 0; i<N-1; i++){
            B[i+1] -= C[i]*A[i+1]/B[i];
            u[i+1] -= u[i]*A[i+1]/B[i];
        }
        u[N-1] *= 1./B[N-1];
        //printf("%f %f %f %f", B[N-1], u[N-2], A[N-2], u[N-1]);
        for (int i = N - 1; i>0; i--){
            u[i-1] = (u[i-1] - u[i]*C[i-1])/B[i-1];
        }
        //printf("%f %f %f \n", B[N-1], C[N-2], u[N-2]);
        break;
    }
};

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
              )
{
    double r, dy, A[M], B[M], C[M];
    int ib2;
    (ib == 0) ? ib2 = 0 : ib2 = ib - 2;
    dy = y[2] - y[1];
    r = x_dcy1/dy;
    A[1] = 0;
    B[1] = x_c1 - x_dcx1/dx - r;
    C[1] = r;
    u[k+m*ib] = x_b1(y[1]) - x_dcx1/dx*u[k+m*ib+ib2];
    for (int i = 2; i<M-2; i++){
        dy = y[i+1] - y[i];
        r = x_dcy1/dy;
        u[k*i+m*ib] = x_b1(y[i]) - x_dcx1/dx*u[k*i+m*ib+ib2];
        A[i] = - r/2;
        B[i] = x_c1 - x_dcx1/dx;
        C[i] = r/2;
    }
    dy = y[M-2] - y[M-3];
    r = x_dcy1/dy;
    A[M-2] = - r;
    B[M-2] = x_c1 - x_dcx1/dx - r;
    C[M-2] = 0;
    u[k*(M-2)+m*ib] = x_b1(y[M-2]) - x_dcx1/dx*u[k*(M-2)+m*ib+ib2];
    for (int i = 1; i<M-2; i++){
        u[k*(i+1)+m*ib] -= u[k*i+m*ib]*A[i+1]/B[i];
        B[i+1] -= C[i]*A[i+1]/B[i];
    }
    u[k*(M-2)+m*ib] *= 1./B[M-2];
    for (int i = M - 3; i>=1; i--){
        u[k*i+m*ib] = (u[k*i+m*ib] - u[k*(i+1)+m*ib]*C[i])/B[i];
    }
}

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
          )
{
    double dx, dy, du, sum_du, a;
    for (int i = 1; i<N-1; i++){
        for (int j = 1; j<M-1; j++){
            u[i*M+j] = 0;
        }
    }
    double A1[N-2], B1[N-2], C1[N-2];
    double A2[M-2], B2[M-2], C2[M-2];
    for (int i = 1; i<N-1; i++){
        u[i*M] = y_b1(x[i]);
        u[i*M+M-1] = y_b2(x[i]);
    }

    dx = x[1] - x[0];
    dy = y[1] - y[0];
    a = 2*(1/dx/dx + 1/dy/dy);

    do {
        go_throw(u, y, x[1] - x[0], 1, M, 0, M, N, x_c1, x_dcx1, x_dcy1, x_b1);
        go_throw(u, y, x[N-1] - x[N-2], 1, M, N-1, M, N, x_c2, x_dcx2, x_dcy2, x_b2);
        go_throw(u, x, y[1] - y[0], N, 1, 0, N, M, y_c1, y_dcy1, y_dcx1, y_b1);
        go_throw(u, x, y[M-1] - y[M-2], N, 1, M-1, N, M, y_c2, y_dcy2, y_dcx2, y_b2);
        sum_du = 0;
        for (int i = 1; i<N-1; i++){
            dx = x[i] - x[i-1];
            for (int j = 1; j<M-1; j++){
                du = omega/a*(1/dx/dx*u[(i+1)*M+j] +
                              1/dx/dx*u[(i-1)*M+j] +
                              1/dy/dy*u[i*M+j+1] +
                              1/dy/dy*u[i*M+j-1] -
                              a*u[i*M+j] -
                              source(x[i], y[j], u[i*M+j])
                              );
                sum_du += fabs(du);
                u[i*M+j] += du;
            }
        }
    printf("%f %f\n", u[N/2*M + M/2], sum_du/N/M/omega);
    } while (sum_du/N/M/omega > eps);
};

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
          )
{
    double dx, dy, du, sum_du = 0, sum_dup, a;
    for (int i = 1; i<N-1; i++){
        for (int j = 1; j<M-1; j++){
            u[i*M+j] = 0;
        }
    }
    u[i1*M + j1] = u1;
    u[i2*M + j2] = u2;
    double A1[N-2], B1[N-2], C1[N-2];
    double A2[M-2], B2[M-2], C2[M-2];
    for (int i = 1; i<N-1; i++){
        u[i*M] = y_b1(x[i]);
        u[i*M+M-1] = y_b2(x[i]);
    }

    dx = x[1] - x[0];
    dy = y[1] - y[0];
    a = 2*(1/dx/dx + 1/dy/dy);

    do {
        for (int j = 0; j<M; j++) u[0*M + j] = (x_b1(y[j]) - x_dcx1*u[1*M + j]/dx)/(x_c1 - x_dcx1/dx);
        for (int j = 0; j<M; j++) u[(N-1)*M + j] = (x_b2(y[j]) + x_dcx2*u[(N-2)*M + j]/dx)/(x_c2 + x_dcx2/dx);
        for (int i = 0; i<N; i++) u[i*M] = (y_b1(x[i]) - y_dcy1*u[i*M + 1]/dx)/(y_c1 - y_dcy1/dy);
        for (int i = 0; i<N; i++) u[i*M + M-1] = (y_b2(y[i]) + y_dcy2*u[i*M + M-2]/dy)/(y_c2 + y_dcy2/dy);
        //go_throw(u, y, x[1] - x[0], 1, M, 0, M, N, x_c1, x_dcx1, x_dcy1, x_b1);
        //go_throw(u, y, x[N-1] - x[N-2], 1, M, N-1, M, N, x_c2, x_dcx2, x_dcy2, x_b2);
        //go_throw(u, x, y[1] - y[0], N, 1, 0, N, M, y_c1, y_dcy1, y_dcx1, y_b1);
        //go_throw(u, x, y[M-1] - y[M-2], N, 1, M-1, N, M, y_c2, y_dcy2, y_dcx2, y_b2);
        sum_dup = sum_du;
        sum_du = 0;
        for (int i = 1; i<N-1; i++){
            dx = x[i] - x[i-1];
            for (int j = 1; j<M-1; j++){
                du = omega/a*(1/dx/dx*u[(i+1)*M+j] +
                              1/dx/dx*u[(i-1)*M+j] +
                              1/dy/dy*u[i*M+j+1] +
                              1/dy/dy*u[i*M+j-1] -
                              a*u[i*M+j]
                              );
                sum_du += fabs(du);
                u[i*M+j] += du;
                if ((i == i1)&&(j == j1)) u[i1*M + j1] = u1;
                if ((i == i2)&&(j == j2)) u[i2*M + j2] = u2;
            }
        }
    printf("%f %f\n", u[N/2*M + M/2], fabs(sum_du - sum_dup));
    } while (fabs(sum_du - sum_dup) > eps);
};
