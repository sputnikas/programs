#include "sde.h"

void right_part(double *y, double t, int n, double *a){
    for (int i = 0; i<n; i++){
        a[i] = - 4*y[n-i-1];
    }
}

void sum_vector(int n, double *y1, double a1, double *y2, double a2, double *result){
    for (int j = 0; j<n; j++){
        result[j] = a1*y1[j] + a2*y2[j];
    }
}

void euler(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result){
    double k1[n];
    rp(y, t, n, k1);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + k1[j]*step;
    }
}

void runge2a(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result){
    double k1[n], k2[n];
    rp(y, t, n, k1);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + k1[j]*step;
    }
    rp(result, t+step, n, k2);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + 1./2*(k1[j] + k2[j])*step;
    }
}

void runge2b(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result){
    double k1[n], k2[n];
    rp(y, t, n, k1);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + k1[j]*step/2;
    }
    rp(result, t+step/2, n, k2);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + k2[j]*step;
    }
}

void runge2c(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result){
    double k1[n], k2[n];
    rp(y, t, n, k1);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + 2./3*k1[j]*step;
    }
    rp(result, t+2*step/3, n, k2);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + 1./4*(k1[j] + 3*k2[j])*step;
    }
}

void runge3(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result){
    double k1[n], k2[n], k3[n];
    rp(y, t, n, k1);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + k1[j]*step/2;
    }
    rp(result, t+step/2, n, k2);
    for (int j = 0; j<n; j++){
        result[j] = y[j] - (k1[j] - 2*k2[j])*step;
    }
    rp(result, t+step, n, k3);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + 1./6*(k1[j] + 4*k2[j] + k3[j])*step;
    }
}

void runge4a(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result){
    double k1[n], k2[n], k3[n], k4[n];
    rp(y, t, n, k1);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + k1[j]*step/2;
    }
    rp(result, t+step/2, n, k2);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + k2[j]*step/2;
    }
    rp(result, t+step/2, n, k3);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + k3[j]*step;
    }
    rp(result, t+step, n, k4);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + 1./6*(k1[j] + 2*k2[j] + 2*k3[j] + k4[j])*step;
    }
}

void runge4b(double *y, int n, double t, double step, void rp(double*, double, int, double*), double *result){
    double k1[n], k2[n], k3[n], k4[n];
    rp(y, t, n, k1);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + k1[j]*step/3;
    }
    rp(result, t+step/3, n, k2);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + (k2[j] - k1[j]/3)*step;
    }
    rp(result, t+2*step/3, n, k3);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + (k1[j] - k2[j] + k3[j])*step;
    }
    rp(result, t+step, n, k4);
    for (int j = 0; j<n; j++){
        result[j] = y[j] + 1./8*(k1[j] + 3*k2[j] + 3*k3[j] + k4[j])*step;
    }
}

//Массив y[n*N]
//Массив y0[n] -- начальные условия в t0
//Решается уравнение y' = rp(y,t) на промежутке t0
void sde(double *y, double *y0, double t0, double tmax, int n, int N, void rp(double*, double, int, double*), int method){
    double step = (tmax-t0)/(N-1);
    double a[n], a_p[n], a_pp[n], a_ppp[n], y1[n], a_n[n];
    double t = t0;
    for (int i = 0; i<n; i++){
        y[i] = y0[i];
    }
    switch (method){
    case 1: // Метод Эйлера
        for (int i = 0; i<(N-1); i++){
            euler(&y[i*n], n, t, step, rp, &y[(i+1)*n]);
            t += step;
        }
        break;
    case 2: // Метод Рунге-Кутты 2-го порядка
        for (int i = 0; i<(N-1); i++){
            runge2a(&y[i*n], n, t, step, rp, &y[(i+1)*n]);
            t += step;
        }
        break;
    case 3: // Метод Рунге-Кутты 2-го порядка
        for (int i = 0; i<(N-1); i++){
            runge2b(&y[i*n], n, t, step, rp, &y[(i+1)*n]);
            t += step;
        }
        break;
    case 4: // Метод Рунге-Кутты 2-го порядка
        for (int i = 0; i<(N-1); i++){
            runge2c(&y[i*n], n, t, step, rp, &y[(i+1)*n]);
            t += step;
        }
        break;
    case 5: // Метод Рунге-Кутты 3-го порядка
        for (int i = 0; i<(N-1); i++){
            runge3(&y[i*n], n, t, step, rp, &y[(i+1)*n]);
            t += step;
        }
        break;
    case 6: // Метод Рунге-Кутты 4-го порядка
        for (int i = 0; i<(N-1); i++){
            runge4a(&y[i*n], n, t, step, rp, &y[(i+1)*n]);
            t += step;
        }
        break;
    case 7: // Метод Рунге-Кутты 4-го порядка
        for (int i = 0; i<(N-1); i++){
            runge4b(&y[i*n], n, t, step, rp, &y[(i+1)*n]);
            t += step;
        }
        break;
    case 8: //Метод Адамса 2-го порядка
        runge4a(&y[0*n], n, t, step, rp, &y[1*n]);
        rp(&y[0*n], t, n, a_p);
        t += step;
        for (int i = 1; i<(N-1); i++){
            rp(&y[i*n], t, n, a);
            for (int j = 0; j<n; j++){
                y1[j] = y[i*n+j] + step*(3*a[j] - a_p[j])/2;
            }
            t+=step;
            rp(y1, t, n, a_n);
            for (int j = 0; j<n; j++){
                y[(i+1)*n+j] = y[i*n+j] + step*(a_n[j] + a[j])/2;
            }
            memcpy(a_p, a, n*sizeof(double));
        }
        break;
    case 9: //Метод Адамса 3-го порядка
        runge4a(&y[0*n], n, t, step, rp, &y[1*n]);
        rp(&y[0*n], t, n, a_p);
        t += step;
        runge4a(&y[1*n], n, t, step, rp, &y[2*n]);
        rp(&y[1*n], t, n, a_pp);
        t += step;
        for (int i = 2; i<(N-1); i++){
            rp(&y[i*n], t, n, a);
            for (int j = 0; j<n; j++){
                y1[j] = y[i*n+j] + step*(23*a[j] - 16*a_p[j] + 5*a_pp[j])/12;
            }
            t+=step;
            rp(y1, t, n, a_n);
            for (int j = 0; j<n; j++){
                y[(i+1)*n+j] = y[i*n+j] + step*(5*a_n[j] + 8*a[j] - a_p[j])/12;
            }
            memcpy(a_pp, a_p, n*sizeof(double));
            memcpy(a_p, a, n*sizeof(double));
        }
        break;
    case 10: //Метод Адамса 4-го порядка
        runge4a(&y[0*n], n, t, step, rp, &y[1*n]);
        rp(&y[0*n], t, n, a_p);
        t += step;
        runge4a(&y[1*n], n, t, step, rp, &y[2*n]);
        rp(&y[1*n], t, n, a_pp);
        t += step;
        runge4a(&y[2*n], n, t, step, rp, &y[3*n]);
        rp(&y[2*n], t, n, a_ppp);
        t += step;
        for (int i = 3; i<(N-1); i++){
            rp(&y[i*n], t, n, a);
            for (int j = 0; j<n; j++){
                y1[j] = y[i*n+j] + step*(55*a[j] - 59*a_p[j] + 37*a_pp[j] - 9*a_ppp[j])/24;
            }
            t+=step;
            rp(y1, t, n, a_n);
            for (int j = 0; j<n; j++){
                y[(i+1)*n+j] = y[i*n+j] + step*(9*a_n[j] + 19*a[j] - 5*a_p[j] + a_pp[j])/24;
            }
            memcpy(a_ppp, a_pp, n*sizeof(double));
            memcpy(a_pp, a_p, n*sizeof(double));
            memcpy(a_p, a, n*sizeof(double));
        }
        break;
    }
}

double norm(double *y, double *y_p, int n){
    double result1 = 0.0;//, result2 = fabs(y[0]);
    for (int i = 0; i<n; i++){
        result1 += fabs(y[i]-y_p[i]);
    }
    return result1;
}

void massive_copy(double *y, double *y_more, int N, int n, int i){
    //y_more[N*k*n];
    //i<k
    //y[N*n]
    //каждые i-е n элементов копируются в y
    for (int j = 0; j<N; j++){
        for (int p = 0; p<n; p++){
            y[j*n+p] = y_more[j*i*n+p];
        }
    }
}

void sde_we(double *y, double *y0, double t0, double tmax, int n, int N, void rp(double*, double, int, double*), int method, double epsilon){
    double *y_c, y_f[n];
    int Nk = N, k = 1;
    y_c = new double[Nk];
    if (method < 20){
        do {
            delete [] y_c;
            y_c = new double[2*Nk*n];
            sde(y_c, y0, t0, tmax, n, Nk, rp, method);
            memcpy(y_f, &y_c[(Nk-1)*n], n*sizeof(double));
            Nk *= 2;
            k *= 2;
            sde(y_c, y0, t0, tmax, n, Nk, rp, method);
        } while (norm(y_f, &y_c[(Nk-1)*n],n)>epsilon);
        massive_copy(y, y_c, N, n, k);
    } else {
        double step = (tmax-t0)/(N-1);
        double a[n], a_p[n], a_pp[n], a_ppp[n], y1[n], a_n[n], d1 = 10000, d2;
        double t = t0;
        for (int i = 0; i<n; i++){
            y[i] = y0[i];
        }
        switch (method){
        case ADAMSSE:
            runge4a(&y[1*n], n, t, step, rp, &y[2*n]);
            rp(&y[0*n], t, n, a_p);
            t += step;
            for (int i = 1; i<(N-1); i++){
                rp(&y[i*n], t, n, a);
                for (int j = 0; j<n; j++){
                    y1[j] = y[i*n+j] + step*(3*a[j] - a_p[j])/2;
                }
                t+=step;
                rp(y1, t, n, a_n);
                for (int j = 0; j<n; j++){
                    y1[j] = y[i*n+j] + step*(a_n[j] + a[j])/2;
                }
                do {
                    memcpy(y_f, y1, n*sizeof(double));
                    rp(y1, t, n, a_n);
                    for (int j = 0; j<n; j++){
                        y1[j] = y[i*n+j] + step*(a_n[j] + a[j])/2;
                    }
                    d2 = d1;
                    d1 = norm(y_f, y1, n);
                } while ((d1>epsilon)&&(d2>d1));
                memcpy(&y[(i+1)*n], y1, n*sizeof(double));
                memcpy(a_p, a, n*sizeof(double));
            }
            break;
        case ADAMSTE:
            runge4a(&y[1*n], n, t, step, rp, &y[2*n]);
            rp(&y[0*n], t, n, a_p);
            t += step;
            runge4a(&y[2*n], n, t, step, rp, &y[3*n]);
            rp(&y[1*n], t, n, a_pp);
            t += step;
            for (int i = 2; i<(N-1); i++){
                rp(&y[i*n], t, n, a);
                for (int j = 0; j<n; j++){
                    y1[j] = y[i*n+j] + step*(23*a[j] - 16*a_p[j] + 5*a_pp[j])/12;
                }
                t+=step;
                rp(y1, t, n, a_n);
                for (int j = 0; j<n; j++){
                    y[(i+1)*n+j] = y[i*n+j] + step*(5*a_n[j] + 8*a[j] - a_p[j])/12;
                }
                do {
                    memcpy(y_f, y1, n*sizeof(double));
                    rp(y1, t, n, a_n);
                    for (int j = 0; j<n; j++){
                        y[(i+1)*n+j] = y[i*n+j] + step*(5*a_n[j] + 8*a[j] - a_p[j])/12;
                    }
                    d2 = d1;
                    d1 = norm(y_f, y1, n);
                } while ((d1>epsilon)&&(d2>d1));
                memcpy(&y[(i+1)*n], y1, n*sizeof(double));
                memcpy(a_pp, a_p, n*sizeof(double));
                memcpy(a_p, a, n*sizeof(double));
            }
            break;
        case ADAMSFE:
            runge4a(&y[1*n], n, t, step, rp, &y[2*n]);
            rp(&y[0*n], t, n, a_p);
            t += step;
            runge4a(&y[2*n], n, t, step, rp, &y[3*n]);
            rp(&y[1*n], t, n, a_pp);
            t += step;
            runge4a(&y[3*n], n, t, step, rp, &y[4*n]);
            rp(&y[2*n], t, n, a_ppp);
            t += step;
            for (int i = 4; i<(N-1); i++){
                rp(&y[i*n], t, n, a);
                for (int j = 0; j<n; j++){
                    y1[j] = y[i*n+j] + step*(55*a[j] - 59*a_p[j] + 37*a_pp[j] - 9*a_ppp[j])/24;
                }
                t+=step;
                rp(y1, t, n, a_n);
                for (int j = 0; j<n; j++){
                    y[(i+1)*n+j] = y[i*n+j] + step*(9*a_n[j] + 19*a[j] - 5*a_p[j] + a_pp[j])/24;
                }
                do {
                    memcpy(y_f, y1, n*sizeof(double));
                    rp(y1, t, n, a_n);
                    for (int j = 0; j<n; j++){
                        y[(i+1)*n+j] = y[i*n+j] + step*(9*a_n[j] + 19*a[j] - 5*a_p[j] + a_pp[j])/24;
                    }
                    d2 = d1;
                    d1 = norm(y_f, y1, n);
                } while ((d1>epsilon)&&(d2>d1));
                memcpy(&y[(i+1)*n], y1, n*sizeof(double));
                memcpy(a_ppp, a_pp, n*sizeof(double));
                memcpy(a_pp, a_p, n*sizeof(double));
                memcpy(a_p, a, n*sizeof(double));
            }
            break;
        }
    }

}
