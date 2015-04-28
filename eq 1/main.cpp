#include <math.h>
#include <stdio.h>

double in_cond(double x){
    return sin(M_PI*x);
}

double v_phase(double u, double x, double t){
    //return sin(u);
    //return u;
    return sqrt(1 + u*u);
    //return sin(u*x);
}

double o_force(double u, double x, double t){
    return -(u - 0.5);
}

int main(){
    double T = 10.;
    double L = 1.;
    int N = 1000;
    int M = 100;
    int k = 20;
    double dt = T/(N-1);
    double dx = L/(M-1);
    double u0[M];
    double umax, umin;
    double xp[M];
    double k1x, k2x, k3x, k4x;
    double k1u, k2u, k3u, k4u;
    
    for (int i = 0; i<M; i++){
        xp[i] = dx*i;
        u0[i] = in_cond(xp[i]);
        if (i==0) { 
            umax = u0[i]; 
            umin = u0[i];
        }
        umax = (u0[i]>umax) ? u0[i] : umax;
        umin = (u0[i]<umin) ? u0[i] : umin;
    }
    
    FILE *gp;
    gp = fopen("eq1.gp", "w");
    fprintf(gp, "set term gif ");
    fprintf(gp, "animate ");
    //fprintf(gp, "optimize ");
    fprintf(gp, "delay 1 ");
    fprintf(gp, "size 600, 600 ");
    fprintf(gp, "background \"#ffffff\" ");
    //fprintf(gp, "crop ");
    fprintf(gp, "\n");
    fprintf(gp, "set output '5.gif' \n ");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set xrange[%e:%e]\n", 0, T);
    fprintf(gp, "set yrange[%e:%e]\n", umin, umax);
    fprintf(gp, "set grid x y\n");
	for (int j = 0; j<N; j++){
        if (j%k == 0) {
            fprintf(gp, "plot '-' using 1:2 with lines lw 2  \n");
            for (int i = 0; i<M; i++) fprintf(gp, "%f %f\n", xp[i], u0[i]);
            fprintf(gp, "end\n");
        }
        for (int i = 0; i<M; i++){
            k1x = v_phase(u0[i], xp[i], j*dt);
            k1u = o_force(u0[i], xp[i], j*dt);
            k2x = v_phase(u0[i] + k1u*dt/2, xp[i] + k1x*dt/2, j*dt + dt/2);
            k2u = o_force(u0[i] + k1u*dt/2, xp[i] + k1x*dt/2, j*dt + dt/2);
            k3x = v_phase(u0[i] + k2u*dt/2, xp[i] + k2x*dt/2, j*dt + dt/2);
            k3u = o_force(u0[i] + k2u*dt/2, xp[i] + k2x*dt/2, j*dt + dt/2);
            k4x = v_phase(u0[i] + k3u*dt, xp[i] + k3x*dt, j*dt + dt);
            k4u = o_force(u0[i] + k3u*dt, xp[i] + k3x*dt, j*dt + dt);
            xp[i] += dt/6*(k1x + 2*k2x + 2*k3x+ k4x);
            u0[i] += dt/6*(k1u + 2*k2u + 2*k3u+ k4u);
        }
	}
	fclose(gp);
    gp = popen("eq1.gp", "r");
    pclose(gp);
}