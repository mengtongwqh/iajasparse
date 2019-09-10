#include <iaja/linalg_matrix.h>

int main() {

    /*  test matrix 1 */
    // unsigned int ia[] = {0, 3, 5, 7};
    // unsigned int ja[] = {0, 1, 2, 0, 1, 0, 2};
    // double a[] = {6.0, 8.0, 8.0, 2.0, 9.0, 5.0, 5.0};

    unsigned int ia[] = {0, 3, 5, 6};
    unsigned int ja[] = {0, 1, 3, 0, 2, 1};
    double a[] = {1, 2, 4, 3, 5, 1};

    unsigned int N = sizeof(ia)/sizeof(unsigned int) - 1;
    unsigned int NNZ = sizeof(a)/sizeof(double);

    iaja::SparseMatrix<double> A(N, 4, NNZ, ia, ja, a);
    iaja::SparseMatrix<double> B = A.transpose();
    std::cout << "Original Matrix: \n";
    std::cout << A << std::endl;
    std::cout << "Transposed Matrix: \n";
    std::cout << B << std::endl;

    return EXIT_SUCCESS;
}
