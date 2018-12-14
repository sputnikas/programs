#pragma once

//#define NONMACROS_TYPE
#ifdef NONMACROS_TYPE
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
