#include <iaja/tests.h>
#include <iaja/incomplete_factor.h>
#include <iaja/iterative_solver.h>

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>


int main(int argc, char **args) {

    // parameters for solver
    if ( argc != 2 )
        std::cerr << "Need 1 argument for grid numbers per each dim." << std::endl;

    std::size_t ngrid = atoi(args[1]);
    unsigned int max_iter = 10000;
    double tol = 1e-10;
    unsigned int max_level_of_fill = 0;

    // construct testing object
    iaja::ImgDenoiseTest testobj(ngrid);

    // build lhs and rhs
    testobj.img_rhs();
    testobj.img_lhs();
    testobj.set_diag_dominant_Mmatrix();

    iaja::SparseIChol* ifac = new iaja::SparseIChol(testobj.get_lhs());
    iaja::PCG solver(ifac, max_iter, tol);

    // test PCG with modified lhs and rhs
    std::ofstream ofm("img.m", std::ofstream::out);

    solver.symbolic_factor(max_level_of_fill);
    std::cout << testobj.x << std::endl;
    if (solver.iterative_solve(testobj.get_rhs(), testobj.x) == EXIT_SUCCESS)
        std::cout << "iterative solver completed: " << solver.get_residual_norm()
            << " after iteration count " << solver.get_iter_count() << std::endl;

    ofm << testobj.x << std::endl;

    // close logging file
    ofm.close();
    // free memory
    delete ifac;

    return EXIT_SUCCESS;
}
