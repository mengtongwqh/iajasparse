#include <cmath>
#include <fstream>
#include "img_matrix.h"

using iaja::ImgMatrixOrthomin;
using iaja::ImgMatrixPCG;

const double FLOATING_POINT_TOLERANCE = 1e-3;

int main(int argc, char* argv[]) {

    // parameters for solver
    if ( argc != 2 ) {
        std::cerr <<
            "Need 1 argument for grid numbers per each dim."
            << std::endl;
        exit(EXIT_FAILURE);
    }

    std::size_t ngrid = atoi(argv[1]);
    unsigned int max_iter = 10000;
    double tol = 1e-3;
    unsigned int max_level_of_fill = 3;
    unsigned int NORTH = 5;
    std::string reorder_method = "natural";

    /*** CONSTRUCT TEST OBJECT ***/

    // ImgMatrixOrthomin testobj(ngrid, NORTH, max_iter, tol);
    ImgMatrixPCG testobj(ngrid, max_iter, tol);

    testobj.img_rhs();
    testobj.img_lhs();
    testobj.set_diag_dominant_Mmatrix();

    time_t start, after_ilu, after_iter;
    start = clock();
    testobj.solver.symbolic_factor(reorder_method, max_level_of_fill);
    after_ilu = clock();
    std::cout << "ILU symbolic factorization time " << std::scientific << double(after_ilu-start)/double(CLOCKS_PER_SEC) << std::setprecision(5) << std::endl;
    
    std::ofstream spfile("sparsity_pattern.txt", std::ofstream::out);
    std::ofstream loffile("level_of_fill.txt", std::ofstream::out);
    testobj.solver.ilu.print_level_of_fill(loffile);
    testobj.lhs.print_compressed(spfile);
    spfile.close();
    loffile.close();
    return EXIT_SUCCESS;

    if (testobj.solver.iterative_solve(testobj.rhs, testobj.x) == EXIT_SUCCESS)
        std::cout << std::scientific << "iterative solver completed with residual norm " << testobj.solver.get_residual_norm()
            << " after iteration count " << testobj.solver.get_iter_count() << std::endl;
    // std::cout << testobj.solver.ilu  << std::endl;
    after_iter = clock();

    // std::cout << testobj.lhs << std::endl;


    for (std::size_t i = 0; i < testobj.x.length(); ++i) {
        if ( fabs(testobj.x[i] - i - 1) > FLOATING_POINT_TOLERANCE ) {
            std::cerr << std::scientific << "Entry " << (i+1) << " of solution has error "  << testobj.x[i]-i-1 << " is wrong!\n";
            exit(EXIT_FAILURE);
        }
    }
    std::cout << "Andre's standard M-matrix test passed\n";

    std::cout << "Iter solve time " << std::scientific << double(after_iter-after_ilu)/double(CLOCKS_PER_SEC) << std::setprecision(5) << std::endl;
    std::cout << "TOTAL TIME: " << std::fixed << double(after_iter-start)/double(CLOCKS_PER_SEC) << std::setprecision(5) << std::endl;

    return EXIT_SUCCESS;
}
