#include "img_matrix.h"
#include <iostream>
#include <fstream>
#include <string>

using std::cout; using std::endl; using std::ofstream;
using iaja::ImgMatrixOrthomin;

int main(int argc, char **args) {

    if ( argc != 2 )
        std::cerr << "Need 1 argument for grid numbers per each dim." << std::endl;

    // parameters for solver
    unsigned int ngrid = atoi(args[1]);
    unsigned int max_iter = 10000;
    double tol = 1e-2;
    unsigned int max_level_of_fill = 1;
    unsigned int n_orth = 5;

    std::string reorder_method = "natural";


    // construct testing object
    // ImgMatrixPCG testobj(ngrid, max_iter, tol);
    ImgMatrixOrthomin testobj(ngrid, n_orth, max_iter, tol);

    // build lhs and rhs
    testobj.img_rhs();
    testobj.img_lhs();

    // symbolic factorization
    testobj.solver.symbolic_factor(reorder_method, max_level_of_fill);

    /*** OUTER ITERATION ***/
    unsigned int outer_itermax = 10;

    for ( unsigned int k_iter = 0; k_iter < outer_itermax; ++k_iter ) {

        // run pcg iterative solver
        testobj.solver.iterative_solve(testobj.rhs, testobj.x);

        std::cout << "iterative solver completed: " << testobj.solver.get_residual_norm()
            << " after iteration count " << testobj.solver.get_iter_count() << "\n";

        // update LHS with current soln
        testobj.img_lhs();
    }

    ofstream os_x("img.m", ofstream::out);
    os_x << testobj.x << endl;
    os_x.close();

    return EXIT_SUCCESS;
}
