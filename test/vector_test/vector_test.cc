#include <iaja/linalg_vector.h>
using namespace iaja;

int main() {

    double a[] = {1, 2, 3, 4, 5};
    unsigned int na = sizeof(a)/sizeof(double);
    FullVector<double> xa(na, a);
    // std::cout << xa << std::endl;
    FullVector<double> xb(std::move(xa));
    FullVector<double> xc(xb);
    std::cout << xc*xb << std::endl;
    std::cout << xb.norm_l2() << std::endl;
//
    // FullVector<double> xb;
    // xb = xa;
    // double alpha = 4.5;
    // xb.saxpy(alpha, xb, xa);
    // std::cout << xa << std::endl;
    // std::cout << xb << std::endl;
//
    unsigned int N = 10, nnz = 5;
    unsigned int iaf[] = {1, 3, 5, 8, 9};
    double af[] = {2.0, 4.0, 5.0, 10.0, -21.0};
    SparseVector<double> sv(N, nnz, iaf, af);
    SparseVector<double> svv(std::move(sv));
    std::cout << sv << std::endl;
    std::cout << svv << std::endl;

    return EXIT_SUCCESS;
}
