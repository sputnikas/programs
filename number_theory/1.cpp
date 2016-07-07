#include <stdio.h>
#include <math.h>

int power(int a, int n)
{
    return (n == 1) ? a : a*power(a, n - 1);
}

int two_power_in_factorial(int p)
{
    int k = 1;
    int s = 0;
    while (p / power(2, k) != 0)
    {
        s += p / power(2, k);
        k++;
    }
    return s;
}

double zeta(double s, int n = 1400)
{
    return (n > 1) ? 1./ pow(n, s) + zeta(s, n - 1) : 1;
}

int main()
{
    double g = 0.577215664901532860;
    double a = zeta(3), a1 = pow(M_PI, 2)/a/g;
    double b = zeta(5), b1 = pow(M_PI, 4)/b/g;
    double c = zeta(7), c1 = pow(M_PI, 6)/c/g;
    printf("%e %e %e\n", a, b, c);
    printf("%e %e %e\n", a1, b1, c1);
    double d = zeta(2), d1 = pow(M_PI, 2)/d;
    double e = zeta(4), e1 = pow(M_PI, 4)/e;
    double f = zeta(6), f1 = pow(M_PI, 6)/f;
    printf("%e %e %e\n", d, e, f);
    printf("%e %e %e\n", d1, e1, f1);
    return 0;
}
