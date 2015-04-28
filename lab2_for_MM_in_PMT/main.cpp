#include "main.h"

char file_directory[80];
const char *name[20] = { "euler",
                         "runge2a",
                         "runge2b",
                         "runge2c",
                         "runge3",
                         "runge4a",
                         "runge4b",
                         "adams2",
                         "adams3",
                         "adams4",
                         "adams2e",
                         "adams3e",
                         "adams4e"
                        };
double par[20];
void (*r_p)(double *, double , int , double *);

void right_part1(double *y, double t, int n, double *a){
    // y' = \sin t
    a[0] = sin(t);
}

void right_part2(double *y, double t, int n, double *a){
    // y'' + 2 \beta y' + \omega_0^2 y = F_0 \sin \omega t
    a[0] = y[1];
    a[1] = - pow(par[0], 2)*y[0] - 2*par[1]*y[1] - par[2]*sin(par[3]*t);
}

void right_part3(double *y, double t, int n, double *a){
    // y'' + 2 \beta y' + (a + b* \cos \omega_0 t) y = F_0 \sin \omega t
    a[0] = y[1];
    a[1] = - (par[0] - 2*par[1]*cos(par[2]*t))*y[0] - 2*par[3]*y[1] - par[4]*sin(par[5]*t);
}

void TestSDE(double *y0, int n, double t0, double tmax, int N, const char *ch){
    int k =10;
    double y[n*N*k], t = t0, step = (tmax-t0)/(N-1);
    sde(y, y0, t0, tmax, n, N, r_p, EULER);
    sde(&y[n*N*1], y0, t0, tmax, n, N, r_p, RUNGESA);
    sde(&y[n*N*2], y0, t0, tmax, n, N, r_p, RUNGESB);
    sde(&y[n*N*3], y0, t0, tmax, n, N, r_p, RUNGESC);
    sde(&y[n*N*4], y0, t0, tmax, n, N, r_p, RUNGET);
    sde(&y[n*N*5], y0, t0, tmax, n, N, r_p, RUNGEFA);
    sde(&y[n*N*6], y0, t0, tmax, n, N, r_p, RUNGEFB);
    sde(&y[n*N*7], y0, t0, tmax, n, N, r_p, ADAMSS);
    sde(&y[n*N*8], y0, t0, tmax, n, N, r_p, ADAMST);
    sde(&y[n*N*9], y0, t0, tmax, n, N, r_p, ADAMSF);
    FILE *f;
    f = fopen(ch, "w");
    fprintf(f, "# ");
    for (int i = 0; i<k; i++){
        fprintf(f, "%s(", name[i]);
        for (int p = 0; p<n; p++){
            fprintf(f, "y%d, ", p);
        }
        fprintf(f, ") ");
    }
    fprintf(f, "\n");
    for (int i = 0; i<N; i++){
        fprintf(f, "%f", t);
        for (int p = 0; p<k; p++){
            for (int j = 0; j<n; j++){
                fprintf(f, " %f", y[N*n*p+i*n+j]);
            }
        }
        t+=step;
        fprintf(f, "\n");
    }
    fclose(f);
}

void TestSDEWE(double *y0, int n, double t0, double tmax, int N, double epsilon, const char *ch){
    int k =13;
    double y[n*N*k], t = t0, step = (tmax-t0)/(N-1);
    sde_we(y, y0, t0, tmax, n, N, r_p, EULER, epsilon);
    sde_we(&y[n*N*1], y0, t0, tmax, n, N, r_p, RUNGESA, epsilon);
    sde_we(&y[n*N*2], y0, t0, tmax, n, N, r_p, RUNGESB, epsilon);
    sde_we(&y[n*N*3], y0, t0, tmax, n, N, r_p, RUNGESC, epsilon);
    sde_we(&y[n*N*4], y0, t0, tmax, n, N, r_p, RUNGET, epsilon);
    sde_we(&y[n*N*5], y0, t0, tmax, n, N, r_p, RUNGEFA, epsilon);
    sde_we(&y[n*N*6], y0, t0, tmax, n, N, r_p, RUNGEFB, epsilon);
    sde_we(&y[n*N*7], y0, t0, tmax, n, N, r_p, ADAMSS, epsilon);
    sde_we(&y[n*N*8], y0, t0, tmax, n, N, r_p, ADAMST, epsilon);
    sde_we(&y[n*N*9], y0, t0, tmax, n, N, r_p, ADAMSF, epsilon);
    sde_we(&y[n*N*10], y0, t0, tmax, n, N, r_p, ADAMSSE, epsilon);
    sde_we(&y[n*N*11], y0, t0, tmax, n, N, r_p, ADAMSTE, epsilon);
    sde_we(&y[n*N*12], y0, t0, tmax, n, N, r_p, ADAMSFE, epsilon);
    FILE *f;
    f = fopen(ch, "w");
    fprintf(f, "# ");
    for (int i = 0; i<k; i++){
        fprintf(f, "%s(", name[i]);
        for (int p = 0; p<n; p++){
            fprintf(f, "y%d, ", p);
        }
        fprintf(f, ") ");
    }
    fprintf(f, "\n");
    for (int i = 0; i<N; i++){
        fprintf(f, "%f", t);
        for (int p = 0; p<k; p++){
            for (int j = 0; j<n; j++){
                fprintf(f, " %f", y[N*n*p+i*n+j]);
            }
        }
        t+=step;
        fprintf(f, "\n");
    }
    fclose(f);
}


int main(int argv, char **argc){
    SetConsoleCP(1251);
    SetConsoleOutputCP(1251);
    printf("%d \n", argv);
    int all_right = 1, all = 0, p, n_params = 0;
    char ch[80];
    double t0, tmax, *y0, epsilon;
    int N, n = 0, with_eps = 0;
    for (int i = 1; i<argv; i++){
        //printf("%s \n", argc[i]);
        if (strncmp(argc[i], "-f", 2)==0) {
            all++;
            strcpy(ch, &(argc[i][2]));
        }
        if (i>=argv) break;
        if (strncmp(argc[i], "-limits", 7)==0) {
            all++;
            i++;    if (i<argv) t0 = atof(argc[i]); else all_right = 0;
            i++;    if (i<argv) tmax = atof(argc[i]); else all_right = 0;
        }
        if (i>=argv) break;
        if (strncmp(argc[i], "-points", 7)==0) {
            all++;
            i++;    if (i<argv) N = atoi(argc[i]); else all_right = 0;
        }
        if (i>=argv) break;
        if (strncmp(argc[i], "-eps", 4)==0) {
            with_eps = 1;
            all++;
            i++;    if (i<argv) epsilon = atof(argc[i]); else all_right = 0;
        }
        if (i>=argv) break;
        if (strncmp(argc[i], "-type", 5)==0){
            all++;
            p = atoi(&(argc[i][5]));
            switch (p){
            case 1:
                n = 1;
                n_params = 0;
                r_p = right_part1;
                break;
            case 2:
                n = 2;
                n_params = 4;
                r_p = right_part2;
                break;
            case 3:
                n_params = 6;
                n = 2;
                r_p = right_part3;
                break;
            }
            y0 = new double[n];
            for (int j = 0; j<n; j++){
                i++;
                if (i<argv) y0[j] = atof(argc[i]); else all_right = 0;
            }
            for (int j = 0; j<n_params; j++){
                i++;
                if (i<argv) par[j] = atof(argc[i]); else all_right = 0;
            }
        }
        if (i>=argv) break;
        if (strncmp(argc[i], "-h", 2)==0){
            all_right = 0;
        }
    }
    if (all < 4){
        all_right = 0;
    }
    if (all_right == 0){
        printf("Signature must be next:\n");
        printf("-fFile -limits min max -points N -type1(or 2, or 3)");
        printf(" y y' .. y^{(n-1)} parameters[0] parameters[1] .. parameters[n_params-1] \n");
        printf("For type1 solves equation y' = sin t. No parameters. \n");
        printf("For type2 solves equation y'' + 2 beta y' + omega_0^2 y = F_0 sin omega t.");
        printf("Parameters omega_0 beta F_0 omega. \n");
        printf("For type3 solves equation y'' + 2 beta y' + (a - 2 b cos omega_0 t) y = F_0 sin omega t.");
        printf("Parameters a b omega_0 beta F_0 omega. \n");
        return 1;
    }
    printf("File %s\n", ch);
    printf("t0 = %f tmax = %f\n", t0, tmax);
    printf("Number points between t0 and tmax %d\n", N);
    printf("Dimension of problem %d\n", n);
    printf("It solves equation:\n");
    switch (p){
    case 1:
        printf("y' = sin t");
        break;
    case 2:
        printf("y'' + 2*%f y' + %f^2 y = %f sin %f t", par[1], par[0], par[2], par[3]);
        break;
    case 3:
        printf("y'' + 2*%f y' + (%f - 2*%f cos %f t) y = %f sin %f t", par[3], par[0], par[1], par[2], par[4], par[5]);
        break;
    }
    if (with_eps == 0){
        TestSDE(y0, n, t0, tmax, N, ch);
    };
    if (with_eps == 1){
        TestSDEWE(y0, n, t0, tmax, N, epsilon, ch);
    }

    return 0;
};
