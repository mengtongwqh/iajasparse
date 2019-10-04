#include <iaja/iaja_config.h>
#include <iaja/iterative_solver.h>
#include <iaja/util.h>

#include <cmath>
#include <fstream>
#include <typeinfo>

#ifdef VERBOSE
    #define REPORT_AVG_RES_NORM do {std::cout \
        << typeid(*this).name() \
        << ": average residual norm = " \
        << sqrt(rnorm2)/N << std::endl;} while (0);
#else
    #define REPORT_AVG_RES_NORM
#endif

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

void IterativeSolverIFactor::symbolic_factor(
        unsigned int max_level_of_fill) {
    assert(ifac != nullptr);
    ifac->analyse(max_level_of_fill);
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

    TIMER_BEGIN

    // do numerical factorization
    if (ifac) ifac->factor();

    // allocate temp vectors
    // r : residual vector
    // rt: (r-tilde) residual vector transformed by M
    // p : search vector
    // q : A*p temp vector for iteration
    const unsigned int N = A.nrow();
    FullVector<double> r(N), rt(N), p(N), q(N);

    // init search and residual vector
    residual(b, x, r);

    if (ifac)
        ifac->solve(r, rt);
    else
        rt = r;

    p = rt;

    // init iteration parameters
    iter_count = 0;
    const double bnorm2 = b*b;
    double rnorm2 = r*r;
    double r_dot_rt = r*rt;

    // iteration loop
    while ( rnorm2 > tol*tol*bnorm2 && iter_count < max_iter ) {

        REPORT_AVG_RES_NORM

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
        if (ifac)
            ifac->solve(r, rt);
        else
            rt = r;

        // update r (dot) rt
        double r_dot_rt_new = r*rt;

        // update search direction p
        double beta = r_dot_rt_new/r_dot_rt;
        p.saxpy(beta, p, rt);

        rnorm2 = r*r;
        r_dot_rt = r_dot_rt_new;
    } // iter loop

    residual_norm = sqrt(rnorm2);

    TIMER_END

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

    TIMER_BEGIN

    if (ifac) ifac->factor();

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

        REPORT_AVG_RES_NORM

        if (ifac)
            ifac->solve(r, rt);
        else
            rt = r;

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
        x.saxpy(alpha, p, x);

        if ( iter_count % interval_residual_recompute ) {
            r.saxpy(-alpha, Ap, r);
        } else {
            residual(b, x, r);
        }

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

    TIMER_END

    if (iter_count < max_iter) {
        return EXIT_SUCCESS;
    } else {
        std::cerr << "Orthomin failed to converge after max_iter " << max_iter << std::endl;
        exit(EXIT_FAILURE);
    }

}

IAJA_NAMESPACE_CLOSE
