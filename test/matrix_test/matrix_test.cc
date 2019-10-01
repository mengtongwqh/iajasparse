#include <cstddef>
#include <iaja/linalg_vector.h>
#include <iaja/linalg_matrix.h>

int main() {

    /*  test matrix 1 */
    // unsigned int ia[] = {0, 3, 5, 7};
    // unsigned int ja[] = {0, 1, 2, 0, 1, 0, 2};
    // double a[] = {6.0, 8.0, 8.0, 2.0, 9.0, 5.0, 5.0};

    /* std::size_t ia[] = {0, 3, 5, 6}; */
    // std::size_t ja[] = {0, 1, 3, 0, 2, 1};
    // double a[] = {1, 2, 4, 3, 5, 1};
/*  */

    std::size_t ia[] = {0, 4, 7, 9, 10, 13};
    std::size_t ja[] = {0, 1, 2, 4, 0, 1, 3, 0, 2, 3, 0, 1, 4};
    double a[] = {1.0, 2.0, 3.0, 5.0, 2.0, 10.0, 1.0, 3.0, 6.0, 5.0, 5.0, 1.0, 4.0};
    std::size_t N = sizeof(ia)/sizeof(std::size_t) - 1;
    std::size_t NNZ = sizeof(a)/sizeof(double);
    iaja::SparseMatrixIaja<double> A(N, 5, NNZ, ia, ja, a);
    
    double b_arr[] = {1, 2, 4, 3, 5};
    iaja::FullVector<double> b(N, b_arr);
    double c_arr[] = {3, 1, 4, 1, 5};
    iaja::FullVector<double> c(N, c_arr);

    std::cout << N << std::endl;

    std::cout << "Original Matrix: \n";
    std::cout << A << std::endl;

    iaja::SparseMatrixIaja<double> B = A.transpose();
    std::cout << "Transposed Matrix: \n";
    std::cout << B << std::endl;

    std::cout << "Operators: \n";
    std::cout << c + A*b << "\n";
    std::cout << c - A*b << "\n";


    return EXIT_SUCCESS;
}
