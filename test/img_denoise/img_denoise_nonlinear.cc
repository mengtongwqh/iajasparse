#include <iaja/tests.h>
#include <iostream>
#include <fstream>
#include <string>

using std::cout; using std::endl; using std::ofstream;

int main(int argc, char **args) {

    if ( argc != 2 )
        std::cerr << "Need 1 argument for grid numbers per each dim." << std::endl;

    // parameters for solver
    std::size_t ngrid = atoi(args[1]);
    unsigned int max_iter = 10000;
    double tol = 1e-10;
    unsigned int max_level_of_fill = 2;
    unsigned int n_orth = 5;


    // construct testing object
    // ImgMatrixPCG testobj(ngrid, max_iter, tol);
    iaja::ImgDenoiseTest testobj(ngrid);

    // build lhs and rhs
    testobj.img_rhs();
    testobj.img_lhs();

    // symbolic factorization
    iaja::IncompleteFactor* ifac = new iaja::SparseIChol(testobj.get_lhs());
    iaja::PCG solver(ifac, max_iter, tol);
    solver.symbolic_factor(max_level_of_fill);

    /*** OUTER ITERATION ***/
    unsigned int outer_itermax = 10;

    for ( unsigned int k_iter = 0; k_iter < outer_itermax; ++k_iter ) {

        // run pcg iterative solver
        solver.iterative_solve(testobj.get_rhs(), testobj.x);

        std::cout << "iterative solver completed: " << solver.get_residual_norm()
            << " after iteration count " << solver.get_iter_count() << "\n";

        // update LHS with current soln
        testobj.img_lhs();
    }

    ofstream os_x("img.m", ofstream::out);
    os_x << testobj.x << endl;
    os_x.close();

    delete ifac;

    return EXIT_SUCCESS;
}
