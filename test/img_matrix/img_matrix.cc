#include "img_matrix.h"
#include <fstream>

IAJA_NAMESPACE_OPEN

ImgMatrixTest::ImgMatrixTest(unsigned int N_dim, unsigned int max_iter, double tol)
    : n(N_dim), N(n*n), lhs(N, N, 5*N), rhs(N), x(N) {
    alpha =
        (n == 16) ? 4e-2 :
        (n == 32) ? 3e-2 :
        (n == 64) ? 1.5e-2 :
        (n == 128) ? 1.2e-2 : 1e-2;
}


void ImgMatrixTest::img_lhs() {

    double h = 1.0 / (n+1);
    double beta = 1e-6;
    lhs.ia[0] = 0;
    unsigned int k = 0;

    for ( unsigned int j = 0; j < n; ++j ) {
        for ( unsigned int i = 0; i < n; ++i ) {

            unsigned int I = j*n + i;
            double a1, a2, a3, a4, a5, a6;

            /* bottom */
            if ( j > 0 ) {
                a1 = (x[I]-x[I-n])*(x[I]-x[I-n]);
                if ( i > 0 ) a1 += (x[I]-x[I-1])*(x[I]-x[I-1]);
                a1 = 1.0/(2*h*sqrt(a1+h*h*beta));
                a2 = (x[I]-x[I-n])*(x[I]-x[I-n]);
                if ( i < n-1 ) a2 += (x[I+1-n]-x[I-n])*(x[I+1-n]-x[I-n]);
                a2 = 1.0/(2*h*sqrt(a2+h*h*beta));
                lhs.ja[k] = I-n; lhs.a[k] = -alpha*(a1+a2); k++;
            }

            /* left */
            if ( i > 0 ) {
                a1 = (x[I]-x[I-1])*(x[I]-x[I-1]);
                if ( j > 0 ) a1 += (x[I]-x[I-n])*(x[I]-x[I-n]);
                a1 = 1.0/(2*h*sqrt(a1+h*h*beta));
                a6 = (x[I]-x[I-1])*(x[I]-x[I-1]);
                if ( j < n-1 ) a6 += (x[I+n-1]-x[I-1])*(x[I+n-1]-x[I-1]);
                a6 = 1.0/(2*h*sqrt(a6+h*h*beta));
                lhs.ja[k] = I-1; lhs.a[k] = -alpha*(a1+a6); k++;
            }

            /* center */
            lhs.ja[k] = I; lhs.a[k] = 0.0; unsigned int kd = k; k++;

            /* right */
            if ( i < n-1 ) {
                a3 = (x[I+1]-x[I])*(x[I+1]-x[I]);
                if ( j > 0 ) a3 += (x[I+1]-x[I+1-n])*(x[I+1]-x[I+1-n]);
                a3 = 1.0/(2*h*sqrt(a3+h*h*beta));
                a4 = (x[I+1]-x[I])*(x[I+1]-x[I]);
                if ( j < n-1 ) a4 += (x[I+n]-x[I])*(x[I+n]-x[I]);
                a4 = 1.0/(2*h*sqrt(a4+h*h*beta));
                lhs.ja[k] = I+1; lhs.a[k] = -alpha*(a3+a4); k++;
            }

            /* top */
            if ( j < n-1 ) {
                a4 = (x[I+n]-x[I])*(x[I+n]-x[I]);
                if ( i < n-1 ) a4 += (x[I+1]-x[I])*(x[I+1]-x[I]);
                a4 = 1.0/(2*h*sqrt(a4+h*h*beta));
                a5 = (x[I+n]-x[I])*(x[I+n]-x[I]);
                if ( i > 0 ) a5 += (x[I+n]-x[I-1+n])*(x[I+n]-x[I-1+n]);
                a5  = 1.0/(2*h*sqrt(a5+h*h*beta));
                lhs.ja[k] = I+n; lhs.a[k] = -alpha*(a4+a5); k++;
            }

            for (unsigned int l = lhs.ia[I]; l < k; l++) {
                if (l != kd) {
                    lhs.a[kd] -= lhs.a[l];
                }
            }
            lhs.a[kd] += 1.0;
            lhs.ia[I+1] = k;
        }
    }
}


void ImgMatrixTest::img_rhs() {

    double h = 1.0/(n+1);

    for ( unsigned int j = 1; j <= n; j++ ) {
        for ( unsigned int i =1 ; i <= n; i++ ) {

            unsigned int I = (j-1)*n+i-1;

            /* exact image value */
            double x = i*h; double y = j*h;
            double r1 = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
            double r2 = sqrt((x-0.5)*(x-0.5)/2 + (y-0.5)*(y-0.5));
            rhs[I] = 0.0;
            if ( r1 < 0.2 ) rhs[I] -= 0.5;
            if ( r2 < 0.3 ) rhs[I] += 1.0;
        }
    }

    /* Add noise to b */
    double b_norm = 0.0;
    unsigned int N = n*n;
    for ( unsigned int i = 0; i < N; i++ ) b_norm += rhs[i]*rhs[i];
    b_norm = sqrt(b_norm);
    srand(42);
    for (unsigned int i = 0; i < N; i++)
        rhs[i] += 2*h*b_norm*(2.0*rand()/RAND_MAX-1);

    /* Copy b to u */
    for (unsigned int i = 0; i < N; i++)
        x[i] = rhs[i];
}


void ImgMatrixTest::set_diag_dominant_Mmatrix() {
/* -------------------------------------------------------------------
 * recondition [a] into a diagonally dominant M-matrix
 * INPUT:
 * N          size of matrix (NxN)
 * ia, ja     sparse pattern of input compressed row matrix
 *
 * INPUT-OUTPUT:
 * a   entries of the original matrix, upon return diagonal will be recomputed
 *
 * OUTPUT:
 * b   rhs of the system Ax = b so as to yield a soln {1,2,3,4,...N}
 * ------------------------------------------------------------------- */

    int isave;

    for ( unsigned int i = 0; i < N; ++i ) {

        isave = -1;
        double temp = 0.0;

        for (unsigned int ii = lhs.ia[i]; ii < lhs.ia[i+1]; ++ii) {
            if ( lhs.ja[ii] == i ) {
                isave = ii;
            } else {
                lhs.a[ii] = -1.0;
                temp += lhs.a[ii];
            }
        }

        if ( isave == -1 ) {
            std::cout << "Error: Row " << i
                << " in ia-ja sparse matrix has a 0 diagonal" << std::endl;
            exit(EXIT_FAILURE);
        }

        lhs.a[isave] = -temp + 0.1;
    }

    for ( unsigned int i = 0; i < N; ++i ) {
        rhs[i] = 0.0;
        for ( unsigned int ii = lhs.ia[i]; ii < lhs.ia[i+1]; ++ii ) {
            rhs[i] += lhs.a[ii] * (lhs.ja[ii] + 1);
        }
    }
}


void ImgMatrixPCG::test_ilu_procedures(unsigned int max_level_of_fill, const std::string& reorder_method) {
    std::cout << "Test ILU factorization with level-of-fill = " << max_level_of_fill
        << " with reordering method \"" << reorder_method << "\"" << std::endl;
    solver.test_ilu_factorization(reorder_method, max_level_of_fill, rhs);
}

IAJA_NAMESPACE_CLOSE
