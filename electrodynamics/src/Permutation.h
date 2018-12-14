#pragma once

#include "Types.h"
#include <iostream>
#include <cstdarg>

void cyclesPermutation(const u32& size, const u32 *a) {
    s32 r = 1;
    char b[size] = {0};
    u32 next;
    for (u32 i = 0; i<size; i++) {
        if (!b[i]) {
            next = i;
            std::cout << "{";
            do {
                if (b[next] == 0) r *= -1;
                next = a[next];
                std::cout << next << ",";
                if (b[next]) {
                    std::cout << std::endl;
                    return;
                }
                b[next] = 1;
            } while (next != i);
            std::cout << "}";
        }
    }
    std::cout << std::endl;
}

void cyclesPermutation(const u32& size, ...) {
    u32 a[size];
    va_list arg;
    va_start(arg, size);
    for (u32 i = 0; i<size; i++) {
        a[i] = va_arg(arg, u32);
    }
    va_end(arg);
    cyclesPermutation(size, a);
}

s32 isOddPermutation(const u32& size, const u32 *a) {
    s32 r = 1;
    bool b[size] = {0};
    u32 next;
    for (u32 i = 0; i<size; i++) {
        if (!b[i]) {
            next = i;
            do {
                if (b[next] == 0) r *= -1;
                next = a[next];
                if (b[next]) return 0;
                b[next] = 1;
            } while (next != i);
        }
    }
    return r;
}

s32 isOddPermutation(const u32& size, ...) {
    u32 a[size];
    va_list arg;
    va_start(arg, size);
    for (u32 i = 0; i<size; i++) {
        a[i] = va_arg(arg, u32);
    }
    va_end(arg);
    return isOddPermutation(size, a);
}

////////////////////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////////////////////

void testPermutation() {
    u32 N = 6;
    u32 a1[N] = {0, 1, 2, 3, 4, 5};
    u32 a2[N] = {0, 1, 2, 3, 5, 4};
    u32 a3[N] = {0, 1, 2, 5, 3, 4};
    u32 a4[N] = {0, 1, 5, 2, 3, 4};
    u32 a5[N] = {0, 1, 5, 3, 2, 4};
    u32 a6[N] = {0, 0, 1, 0, 0, 0};
    cyclesPermutation(N, a6);
    cyclesPermutation(N, a5);
    std::cout << isOddPermutation(N, a1) << std::endl;
    std::cout << isOddPermutation(N, a2) << std::endl;
    std::cout << "list: " << isOddPermutation(N, a3) << std::endl;
    std::cout << "va_list: " << isOddPermutation(N, 0, 1, 2, 5, 3, 4)\
        << std::endl;
    cyclesPermutation(N, a3);
    cyclesPermutation(N, 0, 1, 2, 5, 3, 4);
    std::cout << isOddPermutation(N, a4) << std::endl;
    std::cout << isOddPermutation(N, a5) << std::endl;
    std::cout << isOddPermutation(N, a6) << std::endl;
}


