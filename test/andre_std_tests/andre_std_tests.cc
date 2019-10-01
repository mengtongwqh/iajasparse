#include <iaja/iaja_config.h>
#include <iaja/incomplete_factor.h>
#include <iaja/tests.h>

#include <iostream>
#include <fstream>

using namespace iaja;

struct Parameter {
    /* paramters modified at run time */
    std::string solver_type;
    unsigned int N = 0;
    int max_LoF = -1;
    std::string preconditioner = "ilu";

    /* modify these parameters prior to compile time */
    bool file_output = true;
    unsigned int n_orth = 5;
    unsigned int max_iter = 500;
    const double tol = 1.0e-6;
    const std::string test_method = "aniso";
};


void parse_input(int argc, char *argv[], Parameter& prm) {

    if (argc < 4) {
        std::cerr << "Usage:\n"
            << "(PROG) SolverType GridEachDim MaxLoF (North)"
            << std::endl;
        exit(EXIT_FAILURE);
    }

    prm.preconditioner = argv[1];
    prm.solver_type = argv[2];

    if (prm.solver_type == "pcg" || prm.solver_type == "orthomin") {
        if (argc == 3) {
            prm.N = atoi(argv[3]);
        } else {
            prm.N = atoi(argv[3]);
            prm.max_LoF = atoi(argv[4]);
        }
    } else {
        std::cerr << "Unknown Solver Type" << std::endl;
        exit(EXIT_FAILURE);
    }
}


int main(int argc, char *argv[]) {

    Parameter prm;
    parse_input(argc, argv, prm);

    /* File Output */
    std::ofstream spfile("sparsity_pattern.txt", std::ofstream::out);
    std::ofstream loffile("level_of_fill.txt", std::ofstream::out);
    std::ofstream ilufile("incomplete_factored.txt", std::ofstream::out);
    std::ofstream fullfile("full_pattern.txt", std::ofstream::out);

    EllipticalFDTest testobj(prm.N, prm.test_method);

    IncompleteFactor *ifac;
    if (prm.preconditioner == "ilu")
        ifac = new SparseILU(testobj.get_lhs());
    else if (prm.preconditioner == "ichol")
        ifac = new SparseIChol(testobj.get_lhs());
    else
        exit(EXIT_FAILURE);


    IterativeSolverIFactor *solver;

    if (prm.solver_type == "pcg") {
        prm.max_LoF < 0 ?
            solver = new PCG(testobj.get_lhs(), prm.max_iter, prm.tol):
            solver = new PCG(ifac, prm.max_iter, prm.tol);
    } else if (prm.solver_type == "orthomin") {
        prm.max_LoF < 0 ?
            solver = new Orthomin(testobj.get_lhs(), prm.max_iter, prm.tol, prm.n_orth):
            solver = new Orthomin(ifac, prm.max_iter, prm.tol, prm.n_orth);
    }

    if (prm.max_LoF >= 0) solver->symbolic_factor(prm.max_LoF);
    solver->iterative_solve(testobj.get_rhs(), testobj.x);
    std::cout << solver->get_iter_count() << std::endl;


    /* close output files */
    spfile.close();
    loffile.close();
    ilufile.close();
    fullfile.close();

    /* free memory */
    delete solver;
    delete ifac;

    return EXIT_SUCCESS;
}


