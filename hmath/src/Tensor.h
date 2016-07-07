#ifndef TENSOR_H
#define TENSOR_H

#include "Array.h"

template<typename T> class Tensor : Array<T>
{
private:
    uint32_t rank;
    uint32_t dim;
public:
    Tensor();
    Tensor(const Tensor &a);

    uint32_t get_rank();
    uint32_t get_dim();
    void set_rank(uint32_t r);
    void set_dim(uint32_t d);


}

#endif // TENSOR_H
