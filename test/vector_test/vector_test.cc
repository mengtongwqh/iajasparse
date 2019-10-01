#include <iaja/linalg_vector.h>
using namespace iaja;

int main() {

    // double a[] = {1, 2, 3, 4, 5};
    // FullVector<double>::size_type na = sizeof(a)/sizeof(double);
    // FullVector<double> xa(na, a);
    // std::cout << xa[2] << std::endl;
    // FullVector<double> xb(std::move(xa));
    // FullVector<double> xc(xb);
    // std::cout << xc*xb << std::endl;
    // std::cout << xb.norm_l2() << std::endl;
    // std::cout << xc << std::endl;
    // xc.cumsum();
    // std::cout << "testing iterator...\n";
    // for (FullVector<double>::iterator i = xc.begin(); i != xc.end(); ++i )
        // std::cout << *i << " ";
    // std::cout << std::endl;
    // std::cout << std::endl;
    // for (auto i : xc)
        // std::cout << i << " ";
    // std::cout << std::endl;
    // for (const auto i : xc)
        // std::cout << i << " ";
    // std::cout << std::endl;
    // // std::cout << xc << std::endl;
    // std::cout << "First test completed\n";

    // FullVector<double> xb;
    // xb = xa;
    // double alpha = 4.5;
    // xb.saxpy(alpha, xb, xa);
    // std::cout << xa << std::endl;
    // std::cout << xb << std::endl;

    // unsigned int N = 10, nnz = 5;
    // SparseVector<double>::size_type iaf[] = {1, 3, 5, 8, 9};
    // double af[] = {2.0, 4.0, 5.0, 10.0, -21.0};
    // SparseVector<double> sv(N, nnz, iaf, af);
    // SparseVector<double> svv(std::move(sv));
    // std::cout << sv << std::endl;
    // std::cout << svv << std::endl;

    unsigned int N = 10;

    double *c = new double[N]();
    c[4] = 20; c[5] = 10;
    FullVector<double> cvec(N, std::move(c));
    std::cout << (cvec) << std::endl;
    delete[] c;

    return EXIT_SUCCESS;
}
