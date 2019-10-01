#include <iaja/iaja_config.h>
#include <iaja/linalg_vector.h>
#include <iaja/linalg_matrix.h>
#include <iostream>
using iaja::SizeType;
using iaja::FloatType;
using iaja::FullVector;
using iaja::SparseMatrixIaja;

int main(void) {

    const SizeType N = 6;
    const SizeType NNZ = 20;

    FullVector<SizeType> ia(N+1, new SizeType[N+1] {0, 4, 6, 9, 13, 16, 20});
    FullVector<SizeType> ja(NNZ, new SizeType[NNZ] {0, 2, 3, 5, 1, 3, 0, 2, 4, 0, 1, 3, 5, 2, 4, 5, 0, 3, 4, 5});
    FullVector<FloatType> a (NNZ, new FloatType[NNZ]{
        3.0, -1.0, -1.0, -1.0,
        2.0, -1.0,
        -1.0, 3.0, -1.0,
        -1.0, -1.0, 2.0, -1.0,
        -1.0, 3.0, -1.0,
        -1.0, -1.0, -1.0, 4.0
    });

    SparseMatrixIaja<FloatType> A(N, std::move(ia), std::move(ja), std::move(a));
    FullVector<FloatType> b(N, new FloatType[N]{1, 2, 3, 4, 5, 6});

    std::cout << A << std::endl;
    std::cout << b << std::endl;

    FullVector<FloatType> c;
    c = A*(b+b);
    std::cout << c << std::endl;

    return 0;
}
