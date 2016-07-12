#pragma once

#include <vector>

template <typename T>
struct Element {
    T data;
    std::vector<&Element<T>> to;
};

template <typename T>
struct CompareTree {
    std::vector<Symbol> asbuk;
    Element
};
