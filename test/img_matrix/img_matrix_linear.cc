#include "img_matrix.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>

using iaja::ImgMatrixPCG;

int main(int argc, char **args) {

    // parameters for solver
    if ( argc != 2 )
        std::cerr << "Need 1 argument for grid numbers per each dim." << std::endl;

    std::size_t ngrid = atoi(args[1]);
    unsigned int max_iter = 10000;
    double tol = 1e-10;
    unsigned int max_level_of_fill = 0;
    std::string reorder_method = "natural";


    // construct testing object
    // ImgMatrixPCG testobj(ngrid, max_iter, tol);
    ImgMatrixPCG testobj(ngrid, max_iter, tol);

    // build lhs and rhs
    testobj.img_rhs();
    testobj.img_lhs();
    testobj.set_diag_dominant_Mmatrix();
    // run test for incomplete factorization
    testobj.test_ilu_procedures(max_level_of_fill, reorder_method);

    // test PCG with modified lhs and rhs
    std::ofstream ofs("PCGtest_ilu_pcg_linear.txt", std::ofstream::out);
    std::ofstream ofm("img.m", std::ofstream::out);

    ofs << "original matrix is:" << std::endl;
    ofs << testobj.lhs << std::endl;
    ofs << "RHS is:" << std::endl;
    ofs << testobj.rhs << std::endl;
    ofs << "ILU is:" << std::endl;
    testobj.solver.get_ilu_pattern(ofs);

    testobj.solver.symbolic_factor(reorder_method, max_level_of_fill);
    std::cout << testobj.x << std::endl;
    if (testobj.solver.iterative_solve(testobj.rhs, testobj.x) == EXIT_SUCCESS)
        std::cout << "iterative solver completed: " << testobj.solver.get_residual_norm()
            << " after iteration count " << testobj.solver.get_iter_count() << std::endl;

    ofs << "The soln is:" << std::endl;
    ofs << testobj.x << std::endl;
    ofm << testobj.x << std::endl;
    std::cout << "The soln is: "  << testobj.x <<  std::endl;

    // close logging file
    ofm.close();
    ofs.close();

    return EXIT_SUCCESS;
}
