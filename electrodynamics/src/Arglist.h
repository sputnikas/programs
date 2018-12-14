#pragma once

#define ARG_LIST1(T, a)( T a1 )
#define ARG_LIST2(T, a)( T a2, ARG_LIST1(T, a))
#define ARG_LIST3(T, a)( T a3, ARG_LIST2(T, a))
#define ARG_LIST4(T, a)( T a4, ARG_LIST3(T, a))
#define ARG_LIST5(T, a)( T a5, ARG_LIST4(T, a))
#define ARG_LIST6(T, a)( T a6, ARG_LIST5(T, a))
#define ARG_LIST7(T, a)( T a7, ARG_LIST6(T, a))
#define ARG_LIST8(T, a)( T a8, ARG_LIST7(T, a))
#define ARG_LIST9(T, a)( T a9, ARG_LIST8(T, a))
#define ARG_LIST10(T, a)( T a10, ARG_LIST9(T, a))
#define ARG_LIST11(T, a)( T a11, ARG_LIST10(T, a))

void testArgList1(ARG_LIST4(double, a)) {

}

void testArgList() {
    testArgList1();
}
