
#include <iaja/iaja_config.h>
#include <iaja/incomplete_factor.h>
#include <iaja/iterative_solver.h>
#include <iaja/linalg_matrix.h>
#include <iaja/linalg_vector.h>

#include <iostream>
using std::cout;
using std::endl;
using iaja::SparseMatrixIaja;
using iaja::FullVector;
using iaja::Orthomin;
using iaja::SizeType;
using iaja::FloatType;

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
    const SizeType NNZ = 20;
    const SizeType N = 6;
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
    FloatType* barr = new FloatType[N]{2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

    SparseMatrixIaja<double> A(N, std::move(ia), std::move(ja), std::move(a));
    FullVector<double> b(N, std::move(barr));
    cout << "Orginal matrix is: \n" << A << endl;
    cout << "right hand side matrix is: \n" << b << endl;

    // construct solver object
    unsigned int NORTH = 5;
    iaja::SparseILU ifac(A);
    Orthomin solver(&ifac, 10000, 1e-10, NORTH);
    FullVector<double> x(N);
    solver.symbolic_factor(maxlevel);
    solver.iterative_solve(b, x);
    cout << "The solution is " << endl;
    cout << x << endl;
    cout << "Iteration count is " << solver.get_iter_count()
        << " with iteration residual " << solver.get_residual_norm() << endl;

    return EXIT_SUCCESS;
}
