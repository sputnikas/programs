#include "Integrate.h"

void SwapRows(int i, int j, double *matrix, int n){
    if (i!=j) {
        double temp;
        for (int k = 0; k<n; k++){
            temp = matrix[i*n+k];
            matrix[i*n+k] = matrix[j*n+k];
            matrix[j*n+k] = temp;
        }
    }
};

std::string MatrixToString(int n, double *matrix){
    std::ostringstream str_stream;
    for (int i=0; i<n; i++){
        for (int j=0; j<n; j++){
            str_stream << matrix[i*n+j] << " " ;
        }
        str_stream <<  "\n";
    }
    return str_stream.str();
};

template <class X> std::string MassiveToString(int n, X *matrix){
    std::ostringstream str_stream;
    for (int i=0; i<n; i++){
        str_stream << matrix[i] << " " ;
    }
    return str_stream.str();
};

void LUP(int n, double *matrix, double *l, double *u, int *p){
    for (int i = 0; i<n; i++){
        p[i] = i;
        for (int j = 0; j<n; j++){
            l[i*n+j] = matrix[i*n+j];
            u[i*n+j] = 0;
        }
    }
    double f;
    int temp, kk;
    for (int j = 0; j<n; j++){
        f = l[j*n+j];
        kk = j;
        for (int i = j; i<n; i++){
            if (fabs(l[i*n+j])>fabs(f)) {
                f = l[i*n+j];
                kk = i;
            }
        }
        if (f != 0.0){
            SwapRows(j, kk, l, n);
            if (kk != j){
                temp = p[kk];
                p[kk] = p[j];
                p[j] = temp;
            }
            for (int i = j; i<n; i++){
                for (int k = j+1; k<n; k++){
                    if (i == j) l[k*n+i] = l[k*n+i]/f;
                    if (i>j)   l[k*n+i] = l[k*n+i] - l[j*n+i]*l[k*n+j];
                }
            }
        }
    }
    for (int i = 0; i<n; i++){
        for (int j = i; j<n; j++){
            if (j == i) {
                u[i*n+j] = l[i*n+j];
                l[i*n+j] = 1;
            } else {
                u[i*n+j] = l[i*n+j];
                l[i*n+j] = 0;
            }
        }
    }
};

void SLAE(int n, double *matrix, double *right_part, double *solve){
    double *l, *u, *temp_solve;
    l = new double[n*n];
    u = new double[n*n];
    temp_solve = new double[n];
    int *p;
    p = new int[n];
    LUP(n, matrix, l, u, p);
    double f = 1.0;
    for (int i = 0; i<n; i++) f *= u[i*n+i];
    if (f != 0.0) {
        for (int i = 0; i<n; i++) {
            temp_solve[i] = right_part[p[i]];
            for (int j = 0; j<i; j++) {
                temp_solve[i] -= l[i*n+j]*temp_solve[j];
            }
        }
        for (int i = n-1; i>=0; i--) {
            solve[i] = temp_solve[i];
            for (int j = i + 1; j<n; j++) {
                solve[i] -= u[i*n+j]*solve[j];
            }
            solve[i] /= u[i*n+i];
        }
    } else std::cout << "Матрица вырождена" << std::endl;
    delete [] l;
    delete [] u;
    delete [] p;
}

double func(double x){
    return (x!=0) ? x*sin(1/x) : 0;//(fabs(x)>EPSILON) ? log(x) : 0;//x*x*x;//(fabs(x)>EPSILON) ? log(x) : 0;//log(x);//sqrt(x) : 0; //(x!=0) ? x*sin(1/x) : 0;//x*x*x*x;//x*x;//sin(x);//
}

double Integrate(int n, double *nodes, double *weights){
    double result = 0.0;
    for (int i = 0; i<n; i++){
        result+=weights[i]*func(nodes[i]);
    }
    return result;
};

double IntegrateWithStep(int N, double a, double b, int n, double *weights, double *nodes){
    double result = 0.0;
    double step = (b-a)/N;
    double xmin=a, xmax=xmin+step;
    double *true_nodes, *true_weights;
    true_nodes = new double[n];
    true_weights = new double[n];
    for (int i = 0; i<n; i++){
        true_weights[i] = step*weights[i];
    }
    for (int j = 0; j<N; j++){
        for (int i = 0; i<n; i++){
            true_nodes[i] = ((xmax-xmin)*nodes[i]+xmax+xmin)/2;
        }
        result+=Integrate(n,true_nodes,true_weights);
        xmin+=step;
        xmax+=step;
    }
    return result;
}

double IntegrateLeft(int N, double a, double b){
    double weights[2] = {1, 0};
    double nodes[2] = {-1, 1};
    return IntegrateWithStep(N,a,b,2,weights,nodes);
}

double IntegrateRight(int N, double a, double b){
    double weights[2] = {0, 1};
    double nodes[2] = {-1, 1};
    return IntegrateWithStep(N,a,b,2,weights,nodes);
}

double IntegrateMiddle(int N, double a, double b){
    double weights[2] = {0.5, 0.5};
    double nodes[2] = {-1, 1};
    return IntegrateWithStep(N,a,b,2,weights,nodes);
}

double IntegrateSimpsone(int N, double a, double b){
    double weights[3] = {1./6, 4./6, 1./6};
    double nodes[3] = {-1, 0, 1};
    return IntegrateWithStep(N,a,b,3,weights,nodes);
}

double IntegrateThird(int N, double a, double b){
    double weights[4] = {1./8, 3./8, 3./8, 1./8};
    double nodes[4] = {-1, -1/3, 1/3, 1};
    return IntegrateWithStep(N,a,b,4,weights,nodes);
}

double IntegrateNth(int n, double a, double b){

}

double IntegratePolynomial(int n, double a, double b, double *x){
    double *matrix, *solve;
    matrix = new double[n*n];
    double *right_part;
    right_part = new double[n];
    solve = new double[n];
    for (int i = 0; i<n; i++){
        matrix[i*n] = 1;
        right_part[i] = func(x[i]);
        for (int j = 1; j<n; j++){
            matrix[i*n+j] = matrix[i*n+j-1]*x[i];
        }
    }
    SLAE(n, matrix, right_part, solve);
    //std::cout << MassiveToString(n, solve) << std::endl;
    double min_lim = 1.0, max_lim = 1.0, result = 0.0;
    for (int i = 0; i<n; i++){
        min_lim*=a;
        max_lim*=b;
        result+=solve[i]*(max_lim - min_lim)/(i+1);
    }
    return result;
}

double IntegratePolynomial(int n, double a, double b){
    double x[n], step = (b-a)/(n-1);
    x[0] = a;
    for (int i = 1; i<n; i++){
        x[i] = x[i-1]+step;
    }
    return IntegratePolynomial(n, a, b, x);
}

double IntegrateLagranje(int n, double a, double b, double *x){
    double *a_,*xx;
    a_ = new double[n+1];
    xx = new double[n];
    for (int i = 0; i<n+1; i++) a_[i] = 0.;
    for (int i = 0; i<n; i++){
        xx[i] = 2./(b-a)*(x[i] - (b + a)/2);
    }
    a_[0] = 1.;
    for (int k = 0; k<n; k++){
        for (int i = k+1; i>=0; i--){
            (i==0) ? a_[i] = - a_[i]*xx[k] : a_[i] = a_[i-1] - a_[i]*xx[k];
        }
    }
    //std::cout << MassiveToString(n+1, a_) << std::endl;
    double *b_, *weights, max_lim, min_lim, x_, norma;
    b_ = new double[n];
    weights = new double[n];
    for (int i = 0; i<n; i++){
        b_[0] = 0.;
        max_lim = 1., min_lim = 1., norma = 0., x_ = 1.;
        weights[i] = 0.;
        for (int k = 0; k<n; k++){
            if (fabs(xx[i])>EPSILON) {
                (k!=0) ? b_[k] = (b_[k-1] - a_[k])/xx[i] : b_[k] = - a_[k]/xx[i];
            } else {
                b_[k] = a_[k+1];
            }
            min_lim *= -1.;
            //std::cout << b_[k] << x_ << std::endl;
            norma += b_[k]*x_;
            //std::cout << norma << std::endl;
            x_ *= xx[i];
            weights[i] += b_[k]*(max_lim - min_lim)/(k+1);
        }
        weights[i]/=norma;
        //std::cout << MassiveToString(n, b_) << std::endl;
    }
    //std::cout << "Веса" << MassiveToString(n, weights) << std::endl;
    return (b-a)/2*Integrate(n, x, weights);
}

double IntegrateLagranje(int n, double a, double b){
    double x[n], step = (b-a)/(n-1);
    x[0] = a;
    for (int i = 1; i<n; i++){
        x[i] = x[i-1]+step;
    }
    return IntegrateLagranje(n, a, b, x);
}

int factorial(int n){
    return (n!=0) ? n*factorial(n-1) : 1;
}

int power(int a, int n){
    return (n!=0) ? a*power(a,n-1) : 1;
}

void PolynomLejandreCoef(int n, double *lejandre){
    lejandre[n] = (double) factorial(2*n)/factorial(n)/factorial(n)/power(2, n);
    //std::cout << "Коэффициенты полинома Лежандра     " << lejandre[n] << std::endl;
    int k = 0;
    for (int i = n-1; i>=0; i--){
        if ((n-i)%2 == 0) {
            lejandre[i] = (-1.)*(n-2*k)*(n-2*k-1)/(k+1)/(2*n-2*k-1)/2*lejandre[i+2];
            k++;
        } else {
            lejandre[i] = 0;
        }
    }
    //std::cout << "Коэффициенты полинома Лежандра" << MassiveToString(n+1, lejandre) << std::endl;
}

void PolynomDerivative(int n, double *polynom, double *dpolynom){
    for (int i = 0; i<n; i++){
        dpolynom[i] = polynom[i+1]*(i+1);
    }
    //std::cout << "Производная от полинома " << n << " степени " << MassiveToString(n, dpolynom) << std::endl;
}

double Polynom(int n, double *polynom, double x){
    double result, x_=1.;
    result = polynom[0];
    for (int i = 0; i<n; i++){
        x_*=x;
        result+=polynom[i+1]*x_;
    }
    return result;
}

void PolynomLejandreRoot(int n, double *lejandre, double *root){
    double delta, der[n];
    PolynomDerivative(n, lejandre, der);
    for (int i = 0; i<n; i++){
        root[i] = cos(PI*(4*i+3)/(4*n+2));
        do {
            delta = Polynom(n, lejandre, root[i])/Polynom(n-1, der, root[i]);
            root[i]-=delta;
        } while (fabs(delta)>1e-6);
    }
    //std::cout << Polynom(n-1, der, root[0]) << std::endl;
}

double IntegrateGauss(int n, double a, double b){
    double nodes[n];
    double weights[n];
    double lejandre[n+1];
    PolynomLejandreCoef(n, lejandre);
    PolynomLejandreRoot(n, lejandre, nodes);
    //double d_lejandre[n];
    //double tmp;
    //    PolynomDerivative(n, lejandre, d_lejandre);
    //    for (int j = 0; j<n; j++){
    //        tmp = Polynom(n-1, d_lejandre, nodes[j]);
    //        weights[j] = 2./(1-nodes[j]*nodes[j])/tmp/tmp;
    //    }
    double *matrix;
    matrix = new double[n*n];
    double *right_part;
    right_part = new double[n];
    for (int j = 0; j<n; j++){
        matrix[j] = 1;
        right_part[j] = (j%2==0) ? 2./(j+1) : 0;
        for (int i = 1; i<n; i++){
            matrix[i*n+j] = matrix[i*n+j-n]*nodes[j];
        }
    }
    SLAE(n, matrix, right_part, weights);
    double true_nodes[n];
    for (int i = 0; i<n; i++){
        true_nodes[i] = ((b-a)*nodes[i]+a+b)/2;
    }
    return (b-a)/2*Integrate(n, true_nodes, weights);
}

double IntegrateTanhSinh(int N, double a, double b){
    double infty = 3., step = infty/(N-1), t = - infty, result = 0;
    for (int i = 0; i<N; i++){
        result += (b-a)/2*step*PI*cosh(t)/2/(cosh(PI/2*sinh(t))*cosh(PI/2*sinh(t)))*func((b-a)/2*tanh(PI/2*sinh(t))+(b+a)/2);
        t += step;
    }
    return result;
}

void TestLUP(){
    double a[9] = {2, 7, -6, 8, 2, 1, 7, 4, 2};
    double l[9], u[9];
    int p[3];
    double b[3] = {1, 2, 3};
    double s[3];
    int n = 3;
    LUP(n, a, l, u, p);
    std::cout << MatrixToString(n, a) << std::endl;
    std::cout << MatrixToString(n, l) << std::endl;
    std::cout << MatrixToString(n, u) << std::endl;
    SLAE(n, a, b, s);
    std::cout << MassiveToString(n, s) << std::endl;
    double a1[9] = {1, 2, 0, 3, 4, 4, 5, 6, 3};
    double b1[3] = {3, 7, 8};
    n = 3;
    LUP(n, a1, l, u, p);
    std::cout << MatrixToString(n, a1) << std::endl;
    std::cout << MatrixToString(n, l) << std::endl;
    std::cout << MatrixToString(n, u) << std::endl;
    SLAE(n, a1, b1, s);
    std::cout << MassiveToString(n, s) << std::endl;
    double l1[16],u1[16]; int p1[4];
    double a2[16] = {2, 3, 1, 5 ,6 ,13, 5, 19, 2, 19, 10, 23, 4, 10, 11, 31};
    double b2[4] = {11, 43, 54, 56};
    n = 4;
    LUP(n, a2, l1, u1, p1);
    std::cout << MatrixToString(n, a2) << std::endl;
    std::cout << MatrixToString(n, l1) << std::endl;
    std::cout << MatrixToString(n, u1) << std::endl;
    SLAE(n, a2, b2, s);
    std::cout << MassiveToString(n, s) << std::endl;
    double a3[16] = {1, 1, 1, 1 ,2 ,1, 1, 1, 3, 2, 1, 1, 4, 3, 2, 1};
    double b3[4] = {10, 11, 14, 20};
    n = 4;
    LUP(n, a3, l1, u1, p);
    std::cout << MatrixToString(n, a3) << std::endl;
    std::cout << MatrixToString(n, l1) << std::endl;
    std::cout << MatrixToString(n, u1) << std::endl;
    SLAE(n, a3, b3, s);
    std::cout << MassiveToString(n, s) << std::endl;
}

template <class X> void Output(const char* method, X result, X time){
    std::cout.width(30);
    std::cout << method;
    std::cout.width(15);
    std::cout << result;
    std::cout.width(15);
    std::cout << time;
    std::cout << std::endl;
}

double IntegrateAdaptive(int N, double a, double b, int n, double int_f(int, double, double)){
    double step = (b-a)/(N-1), x = a, result = 0.;
    for (int i = 0; i<N; i++){
        result+=int_f(n, x, x+step);
        x+=step;
    }
    return result;
}

void TestIntegrate(){
    double maxlim = 1;
    int n = 16;
    double result[100], time[100];
    using namespace std::chrono;
    steady_clock::time_point start, finish;
    start = steady_clock::now(); result[0] = IntegrateLeft(n, 0, maxlim);           finish = steady_clock::now(); time[0]  = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[1] = IntegrateRight(n, 0, maxlim);          finish = steady_clock::now(); time[1]  = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[2] =  IntegrateMiddle(n, 0, maxlim);        finish = steady_clock::now(); time[2]  = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[3] =  IntegrateSimpsone(n, 0, maxlim);      finish = steady_clock::now(); time[3]  = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[4] =  IntegrateThird(n, 0, maxlim);         finish = steady_clock::now(); time[4]  = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[5] =  IntegratePolynomial(4, 0, maxlim);    finish = steady_clock::now(); time[5]  = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[6] =  IntegratePolynomial(8, 0, maxlim);    finish = steady_clock::now(); time[6]  = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[7] =  IntegratePolynomial(51, 0, maxlim);   finish = steady_clock::now(); time[7]  = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[8] =  IntegrateLagranje(2, 0, maxlim);      finish = steady_clock::now(); time[8]  = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[9] = IntegrateLagranje(3, 0, maxlim);       finish = steady_clock::now(); time[9]  = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[10] = IntegrateLagranje(4, 0, maxlim);      finish = steady_clock::now(); time[10] = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[11] = IntegrateLagranje(8, 0, maxlim);      finish = steady_clock::now(); time[11] = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[12] = IntegrateLagranje(20, 0, maxlim);     finish = steady_clock::now(); time[12] = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[13] = IntegrateGauss(2, 0, maxlim);         finish = steady_clock::now(); time[13] = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[14] = IntegrateGauss(3, 0, maxlim);         finish = steady_clock::now(); time[14] = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[15] = IntegrateGauss(4, 0, maxlim);         finish = steady_clock::now(); time[15] = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[16] = IntegrateGauss(8, 0, maxlim);         finish = steady_clock::now(); time[16] = (duration_cast<duration<double>>(finish - start)).count();
    start = steady_clock::now(); result[17] = IntegrateGauss(15, 0, maxlim);        finish = steady_clock::now(); time[17] = (duration_cast<duration<double>>(finish - start)).count();
    //start = steady_clock::now(); result[18] = IntegrateTanhSinh(5, 0, maxlim);      finish = steady_clock::now(); time[18] = (duration_cast<duration<double>>(finish - start)).count();
    //start = steady_clock::now(); result[19] = IntegrateTanhSinh(10, 0, maxlim);     finish = steady_clock::now(); time[19] = (duration_cast<duration<double>>(finish - start)).count();
    //start = steady_clock::now(); result[20] = IntegrateTanhSinh(15, 0, maxlim);     finish = steady_clock::now(); time[20] = (duration_cast<duration<double>>(finish - start)).count();
    //start = steady_clock::now(); result[21] = IntegrateTanhSinh(20, 0, maxlim);     finish = steady_clock::now(); time[21] = (duration_cast<duration<double>>(finish - start)).count();
    //start = steady_clock::now(); result[22] = IntegrateTanhSinh(13, 0, maxlim);     finish = steady_clock::now(); time[22] = (duration_cast<duration<double>>(finish - start)).count();
    Output("Метод", "Результат", "Время");
    Output("Левые прямоугольники",          result[0],  time[0]);
    Output("Правые прямоугольники",         result[1],  time[1]);
    Output("Средние прямоугольники",        result[2],  time[2]);
    Output("Метод Симпсона",                result[3],  time[3]);
    Output("Метод 3 порядка",               result[4],  time[4]);
    Output("Полиномом 3 степени",           result[5],  time[5]);
    Output("Полиномом 7 степени",           result[6],  time[6]);
    Output("Полиномом 50 степени",          result[7],  time[7]);
    Output("Полиномом Лагранжа 1 степени",  result[8],  time[8]);
    Output("Полиномом Лагранжа 2 степени",  result[9],  time[9]);
    Output("Полиномом Лагранжа 3 степени",  result[10], time[10]);
    Output("Полиномом Лагранжа 7 степени",  result[11], time[11]);
    Output("Полиномом Лагранжа 19 степени", result[12], time[12]);
    Output("Квадратура Гаусса 2 степени",   result[13], time[13]);
    Output("Квадратура Гаусса 3 степени",   result[14], time[14]);
    Output("Квадратура Гаусса 4 степени",   result[15], time[15]);
    Output("Квадратура Гаусса 8 степени",   result[16], time[16]);
    Output("Квадратура Гаусса 15 степени",  result[17], time[17]);
    //Output("TanhSinh по 5 точкам",          result[18], time[18]);
    //Output("TanhSinh по 10 точкам",         result[19], time[19]);
    //Output("TanhSinh по 15 точкам",         result[20], time[20]);
    //Output("TanhSinh по 20 точкам",         result[21], time[21]);
    //Output("TanhSinh по 13 точкам",         result[22], time[22]);
    //std::cout << "Полиномом Лагранжа 7 степени " << IntegrateLagranje(8, 0, maxlim) << std::endl;
    //std::cout << "Полиномом Лагранжа 50 степени " << IntegrateLagranje(51, 0, maxlim) << std::endl;

}

void TestLejandre(int n){
    double lejandre[n+1];
    double der[n];
    PolynomLejandreCoef(n, lejandre);
//    for (int i = 0; i<n+1; i++){
//        lejandre[i] = 1;
//    }
    PolynomDerivative(n, lejandre, der);
}
