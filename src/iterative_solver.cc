#include <iaja/iaja_config.h>
#include <iaja/iterative_solver.h>
#include <cmath>
#include <fstream>

IAJA_NAMESPACE_OPEN

/* ============================================ *
 *              IterativeSolver                 *
 * ============================================ */
void IterativeSolver::residual(
        const FullVector<double>& b, const FullVector<double>& x,
        FullVector<double>& res) {
    // res = b - A*x
    res.multiply(A, x); // res = A*x
    res.minus(b, res); // res = b - res
}


/* ============================================ *
 *           IterativeSolverIFactor             *
 * ============================================ */

void IterativeSolverIFactor::symbolic_factor(const std::string& reorder_method,
        unsigned int max_level_of_fill) {
    // call ILU reordering and analyse
    ifac->reorder(reorder_method);
    ifac->analyse(max_level_of_fill);
}


void IterativeSolverIFactor::test_incomplete_factorization(const std::string& reorder_method,
        unsigned int max_level_of_fill, const FullVector<double>& rhs) {

    std::ofstream ofs("IterSolvertest_ILU.txt", std::ofstream::out);

    ofs << "original matrix is:" << std::endl;
    ofs << A << std::endl;
    ofs << "RHS is:" << std::endl;
    ofs << rhs << std::endl;

    ofs << "computing reordering with method " << reorder_method << std::endl;
    ifac->reorder(reorder_method);
    ofs << "computing ILU with level " << max_level_of_fill << std::endl;
    ifac->analyse(max_level_of_fill);
    ifac->print_level_of_fill(ofs);
    ofs << "computing numerical factorization" << std::endl;
    ifac->factor();
    ofs << *ifac << std::endl;
    FullVector<double> x(A.nrow());
    ifac->solve(rhs, x);
    ofs << x << std::endl;

    ofs.close();
}

/* ============================================ *
 *                PCG SOLVER                    *
 * ============================================ */

int PCG::iterative_solve(const FullVector<double>& b, FullVector<double>& x) {
/* ---------------------------------------
 * INPUT
 *   b: rhs of the system
 *   x: initial guess
 * OUTPUT
 *   x: converged soln
 * --------------------------------------- */

    assert(A.nrow() == A.ncol());

    // do numerical factorization
    ifac->factor();

    // allocate temp vectors
    // r : residual vector
    // rt: (r-tilde) residual vector transformed by M
    // p : search vector
    // q : A*p temp vector for iteration
    const unsigned int N = A.nrow();
    FullVector<double> r(N), rt(N), p(N), q(N);

    // init search and residual vector
    residual(b, x, r);
    ifac->solve(r, rt);
    p = rt;

    // init iteration parameters
    iter_count = 0;
    const double bnorm2 = b*b;
    double rnorm2 = r*r;
    double r_dot_rt = r*rt;

    // iteration loop
    while ( rnorm2 > tol*tol*bnorm2 && iter_count < max_iter ) {

        std::cout << sqrt(rnorm2)/N << std::endl;

        ++iter_count;

        // search step length along p
        q.multiply(A, p);
        double alpha = r_dot_rt / (p*q);

        // update approximate soln x
        x.saxpy(alpha, p, x);

        // update residual r
        if ( iter_count % interval_residual_recompute ) {
            r.saxpy(-alpha, q, r);
        } else {
            residual(b, x, r);
        }

        // solve for r-tilde
        ifac->solve(r, rt);

        // update r (dot) rt
        double r_dot_rt_new = r*rt;

        // update search direction p
        double beta = r_dot_rt_new/r_dot_rt;
        p.saxpy(beta, p, rt);

        rnorm2 = r*r;
        r_dot_rt = r_dot_rt_new;
    } // iter loop

    residual_norm = sqrt(rnorm2);

    if (iter_count < max_iter) {
        return EXIT_SUCCESS;
    } else {
        std::cerr << "PCG failed to converge after max_iter " << max_iter << std::endl;
        exit(EXIT_FAILURE);
    }
}


/* ============================================ *
 *           ORTHOMIN SOLVER                    *
 * ============================================ */
int Orthomin::iterative_solve(const FullVector<double>& b, FullVector<double>& x) {

    assert(A.nrow() == A.ncol());

    ifac->factor();
    const unsigned int N = A.nrow();

    // allocate temp vectors
    FullVector<double> Ap(N), Art(N), p(N), rt(N), r(N);
    residual(b, x, r);

    // initiate iteration parameters
    iter_count = 0;
    unsigned int k = 0;
    double bnorm2 = b*b;
    double rnorm2 = bnorm2;
    std::vector< FullVector<double> > pk, Apk;
    pk.reserve(k_orth); Apk.reserve(k_orth);

    while ( rnorm2 > bnorm2*tol*tol && iter_count < max_iter ) {

        std::cout << sqrt(rnorm2)/N << std::endl;

        ifac->solve(r, rt);

        p = rt;
        Ap.multiply(A, p);
        Art.multiply(A, rt);


        for (unsigned int i = 0; i < k; ++i) {
            double beta = - (Art*Apk[i])/(Apk[i]*Apk[i]);
            p.saxpy(beta, pk[i], p);
            Ap.saxpy(beta, Apk[i], Ap);
        }

        pk.push_back(p);
        Apk.push_back(Ap);

        double alpha = (r*Ap)/(Ap*Ap);
        x.saxpy(alpha , p, x);
        r.saxpy(-alpha, Ap, r);

        // restart orthomin(k) after k_orth steps
        if ( ++k == k_orth ) {
            pk.clear();
            Apk.clear();
            k = 0;
        }

        rnorm2 = r*r;
        ++iter_count;
    }

    residual_norm = sqrt(rnorm2);
    if (iter_count < max_iter) {
        return EXIT_SUCCESS;
    } else {
        std::cerr << "Orthomin failed to converge after max_iter " << max_iter << std::endl;
        exit(EXIT_FAILURE);
    }
}

IAJA_NAMESPACE_CLOSE
