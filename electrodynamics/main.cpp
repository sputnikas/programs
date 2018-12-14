#include "src/Array.h"
#include "src/Vector.h"
#include "src/MathFunction.h"
#include "src/Permutation.h"
#include <cstdio>
#include <cstdarg>


int prod(int num, const int& n, ...) {
    int r = 1, t;
    va_list argptr;
	va_start(argptr, n);
	for(; num; num--) {
		t = va_arg(argptr, int);
		r = (r*t)%n;
	}
	va_end(argptr);
	return r;
}

void testProd(){
    std::cout << prod(1, 10, 2) << std::endl;
    std::cout << prod(3, 10, 1, 2, 3) << std::endl;
    std::cout << prod(5, 10, 1, 2, 3, 4, 5) << std::endl;
}

int main() {
    //testMathFunction();
    //testProd();
    //testPermutation();
    testVector();
    return 0;
}
