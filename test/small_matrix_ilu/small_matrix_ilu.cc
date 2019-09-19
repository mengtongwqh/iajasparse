
#include <iaja/ilu.h>
#include <iaja/linalg_matrix.h>
#include <iaja/linalg_vector.h>

#include <cstddef>
#include <iostream>
using std::cout;
using std::endl;
using iaja::SparseMatrixIaja;
using iaja::FullVector;
using iaja::SparseILU;

/* const int N does not work for 
 * (int[N]){elmt1, elmt2, ..., elmtN}*/

int main(void) {

    const unsigned int maxlevel = 1;

    /* allocate space */
    /* double *x = (double*)malloc( N*sizeof(double) );  [> soln <] */


    /* set up test matrix */

    /*  test matrix 1 */
    // unsigned int ia[] = {0,3,5,7};
    // unsigned int ja[] = {0,1,2,0,1,0,2};
    // double a[] = {6.0, 8.0, 8.0, 2.0, 9.0, 5.0, 5.0};
    // double b[] = {2.0, 3.0, 4.0};

    /* test matrix 2 */
    std::size_t ia[] = {0, 4, 6, 9, 13, 16, 20};
    std::size_t ja[] = {0, 2, 3, 5, 1, 3, 0, 2, 4, 0, 1, 3, 5, 2, 4, 5, 0, 3, 4, 5};
    double a[] = {
        3.0, -1.0, -1.0, -1.0,
        2.0, -1.0,
        -1.0, 3.0, -1.0,
        -1.0, -1.0, 2.0, -1.0,
        -1.0, 3.0, -1.0,
        -1.0, -1.0, -1.0, 4.0
    };
    double barr[] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

    std::size_t N = sizeof(ia)/sizeof(std::size_t) - 1;
    std::size_t NNZ = sizeof(a)/sizeof(double);
    SparseMatrixIaja<double> A(N, N, NNZ, ia, ja, a);
    FullVector<double> b(N, barr);
    cout << "Orginal matrix is: \n" << A;

    // symbolic factorization
    SparseILU ilu(A);
    ilu.reorder("natural");
    ilu.analyse(maxlevel);
    cout << "\nAfter brute-force symbolic factorization:\n";
    ilu.print_level_of_fill(cout);

    // numerical factorization
    ilu.factor();
    cout << "\nAfter numerical factorization:\n";
    cout << ilu;

    // solve
    FullVector<double> x(N);
    ilu.solve(b, x);
    cout << "\nThe solution is:\n";
    cout << x << "\n";

    // test Ax = b
    FullVector<double> f(N);
    f.multiply(A, x);
    cout << f << endl;

    return EXIT_SUCCESS;
}
