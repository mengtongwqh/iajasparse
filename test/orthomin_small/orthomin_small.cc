
#include <iaja/ilu.h>
#include <iaja/linalg_matrix.h>
#include <iaja/linalg_vector.h>
#include <iaja/iterative_solver.h>

#include <iostream>
using std::cout;
using std::endl;
using iaja::SparseMatrix;
using iaja::FullVector;
using iaja::Orthomin;

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
    SparseMatrix<double> A(N, N, NNZ, ia, ja, a);
    FullVector<double> b(N, barr);
    cout << "Orginal matrix is: \n" << A << endl;
    cout << "right hand side matrix is: \n" << b << endl;

    // construct solver object
    unsigned int NORTH = 2;
    Orthomin solver(A, NORTH, 10000, 1e-10);
    FullVector<double> x(N);
    solver.symbolic_factor("natural", maxlevel);
    solver.iterative_solve(b, x);
    cout << "The solution is " << endl;
    cout << x << endl;
    cout << "Iteration count is " << solver.get_iter_count()
        << " with iteration residual " << solver.get_residual_norm() << endl;

    return EXIT_SUCCESS;
}
