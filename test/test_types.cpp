#include <iostream>
#include <chrono>
#include <typeinfo>

//#define MACROS_TYPE
#ifdef MACROS_TYPE
typedef unsigned char          u8;
typedef unsigned short int     u16;
typedef unsigned long int      u32;
typedef unsigned long long int u64;
typedef signed char            s8;
typedef signed short int       s16;
typedef signed long int        s32;
typedef signed long long int   s64;
typedef float                  f32;
typedef double                 f64;
#else
#define u8   unsigned char
#define u16  unsigned short int
#define u32  unsigned long int
#define u64  unsigned long long int
#define s8   signed char
#define s16  signed short int
#define s32  signed long int
#define s64  signed long long int
#define f32  float
#define f64  double
#endif

template <typename T>
void printsize(){
    std::cout << "sizeof " << typeid(T).name() << " " << sizeof(T) << std::endl;
}

template <typename T>
void printsize(T &a){
    std::cout << "sizeof " << typeid(a).name() << " " << sizeof(a) << std::endl;
}

int main() {
    f64 a = 1;
    std::cout << a << a + 1.1 << a + 2.1 << std::endl;
    bool b[16];
    std::cout << sizeof(b) << std::endl;
    printsize(b);
    printsize<unsigned char>();
    printsize<signed char>();
    printsize<short int>();
    printsize<unsigned short int>();
    printsize<signed short int>();
    printsize<int>();
    printsize<unsigned int>();
    printsize<signed int>();
    printsize<long int>();
    printsize<unsigned long int>();
    printsize<signed long int>();
    printsize<long long int>();
    printsize<unsigned long long int>();
    printsize<signed long long int>();
    printsize<float>();
    printsize<double>();
    printsize<long double>();
    printsize<wchar_t>();
    return 0;
}
