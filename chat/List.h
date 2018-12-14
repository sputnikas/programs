#pragma once

#include <cstdio>

template <typename T>
struct ListElement {
    T data;
    ListElement<T> *to;
    ListElement<T> *from;

    ListElement() : to(0), from(0) {};
    ListElement(const T &a) : data(a), to(0), from(0) {};
};

template <typename T>
struct  List {
private:
    unsigned length;
    ListElement<T> *_first = 0;

    ListElement<T> *_get(unsigned i);
public:
    List() : length(0), _first(0) {};
    List(const T &a);
    List(T *a, const unsigned &n);
    ~List();

    void add(const T &a);
    void add(const T &a, const unsigned &i);
    void add(T *a, const unsigned &n, const unsigned &i);
    void set(const T &a, const unsigned &i);
    T last();
    T first();
    T get(const unsigned &i);
    unsigned size();
    void remove(const unsigned &i);
};

template <typename T>
List<T>::List(const T &a) {
    ListElement<T> *ptr = new ListElement<T>(a);
    _first = ptr;
    length = 1;
}

template <typename T>
List<T>::List(T *a, const unsigned &n) : length(n) {
    ListElement<T> *ptr, *ptr2;
    for (unsigned i = 0; i<n; i++) {
        ptr = new ListElement<T>;
        ptr->data = a[i];
        //printf("%d ", a[i]);
        ptr->from = (i == 0) ? 0 : ptr2;
        if (i > 0) ptr2->to = ptr;
        ptr2 = ptr;
        if (i == 0) _first = ptr;
    }
}

template <typename T>
List<T>::~List() {
    ListElement<T> *ptr, *ptr2 = _first;
    length = 0;
    while (ptr2 != 0) {
        ptr = ptr2;
        ptr2 = ptr->to;
        delete ptr;
    }
}

template <typename T>
ListElement<T> *List<T>::_get(unsigned i) {
    ListElement<T> *ptr = _first;
    while ((i--) && (ptr->to != 0)) {
        ptr = ptr->to;
    }
    return ptr;
}

template <typename T>
T List<T>::first() {
    return _first->data;
}

template <typename T>
T List<T>::last() {
    return _get(length)->data;
}

template <typename T>
T List<T>::get(const unsigned &i) {
    return _get(i)->data;
}

template <typename T>
unsigned List<T>::size() {
    return length;
}

template <typename T>
void List<T>::add(const T &a, const unsigned &i) {
    ListElement<T> *ptr = _get(i), *ptr2;
    length++;
    ptr2 = new ListElement<T>(a);
    if (length > 0) {
        if (i == 0) {
            ptr2->to = ptr;
            ptr->from = ptr2;
            _first = ptr2;
        }
        else {
            ptr2->to = ptr;
            ptr2->from = ptr->from;
            if (ptr->from != 0) ptr->from->to = ptr2;
            ptr->from = ptr2;
        }
    }
    else {
        _first = ptr2;
    }
}

template <typename T>
void List<T>::add(const T &a) {
    ListElement<T> *ptr = _get(length), *ptr2;
    length++;
    ptr2 = new ListElement<T>(a);
    ptr2->from = ptr;
    ptr->to = ptr2;
}

template <typename T>
void List<T>::add(T *a, const unsigned &n, const unsigned &i) {
    for (unsigned j = 0; j < n; j++)
        add(a[j], i+j);
}

template <typename T>
void List<T>::set(const T &a, const unsigned &i) {
    if (i < length) {
        ListElement<T> *ptr = _get(i);
        ptr->data = a;
    }
}

template <typename T>
void List<T>::remove(const unsigned &i) {
    if (i < length) {
        ListElement<T> *ptr = _get(i);
        if (ptr->from != 0) ptr->from->to = ptr->to;
        else _first = ptr->to;
        if (ptr->to != 0) ptr->to->from = ptr->from;
        delete ptr;
        length--;
    }
}

template <typename T>
void toConsole(List<T> a, const char *c = "\n") {
    printf("[ ");
    for (unsigned i = 0; i < a.size(); i++) {
        printf("%d ", a.get(i));
    }
    printf("]%s", c);
}

void testList(){
    int b[6] = {1, 2, 4, 6, 5, 3};
    int c[3] = {7, 8, 9};
    List<int> a;
    a = List<int>(b, 6); toConsole<int>(a);
    a.add(10);           toConsole<int>(a);
    a.add(11, 3);        toConsole<int>(a);
    a.add(12, 0);        toConsole<int>(a);
    a.add(c, 3, 7);      toConsole<int>(a);
    a.remove(2);         toConsole<int>(a);
    a.set(13, 5);        toConsole<int>(a);
    printf("%d %d\n", a.first(), a.last());
}
