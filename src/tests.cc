#include <iaja/tests.h>
#include <fstream>

IAJA_NAMESPACE_OPEN

/* ==================================================== *
 *            IMAGE DENOISING TEST PROBLEM              *
 * ==================================================== */

ImgMatrixTest::ImgMatrixTest(size_type N_dim)
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
    size_type k = 0;

    for ( size_type j = 0; j < n; ++j ) {
        for ( size_type i = 0; i < n; ++i ) {

            size_type I = j*n + i;
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
            lhs.ja[k] = I; lhs.a[k] = 0.0; size_type kd = k; k++;

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

            for (size_type l = lhs.ia[I]; l < k; l++) {
                if (l != kd) {
                    lhs.a[kd] -= lhs.a[l];
                }
            }
            lhs.a[kd] += 1.0;
            lhs.ia[I+1] = k;
        }
    }

    lhs.compress_storage();
}


void ImgMatrixTest::img_rhs() {

    double h = 1.0/(n+1);

    for ( size_type j = 1; j <= n; j++ ) {
        for ( size_type i =1 ; i <= n; i++ ) {

            size_type I = (j-1)*n+i-1;

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
    size_type N = n*n;
    for ( size_type i = 0; i < N; i++ ) b_norm += rhs[i]*rhs[i];
    b_norm = sqrt(b_norm);
    srand(42);
    for (size_type i = 0; i < N; i++)
        rhs[i] += 2*h*b_norm*(2.0*rand()/RAND_MAX-1);

    /* Copy b to u */
    for (size_type i = 0; i < N; i++)
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

    for ( size_type i = 0; i < N; ++i ) {

        isave = -1;
        double temp = 0.0;

        for (size_type ii = lhs.ia[i]; ii < lhs.ia[i+1]; ++ii) {
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

    for ( size_type i = 0; i < N; ++i ) {
        rhs[i] = 0.0;
        for ( size_type ii = lhs.ia[i]; ii < lhs.ia[i+1]; ++ii ) {
            rhs[i] += lhs.a[ii] * (lhs.ja[ii] + 1);
        }
    }
}


void ImgMatrixPCG::test_ilu_procedures(unsigned int max_level_of_fill, const std::string& reorder_method) {
    std::cout << "Test ILU factorization with level-of-fill = " << max_level_of_fill
        << " with reordering method \"" << reorder_method << "\"" << std::endl;
    solver.test_ilu_factorization(reorder_method, max_level_of_fill, rhs);
}



/* ==================================================== *
 *           3D DIFFUSION TEST PROBLEM                  *
 * ==================================================== */

FDGrid::FDGrid(size_type n_dim, const std::string& method):
    nx(n_dim), ny(n_dim), nz(n_dim), n(nx*ny*nz), rhs(n, 0.0), x(n, 0.0), lhs(), idiag(n, 0.0) {

    if (method == "ascend") {
        set_ascend();
    } else if (method == "aniso") {
        set_anisotropic_K();
    } else {
        std::cerr << "Unknown method of constructing LHS in FDGrid\n";
        exit(EXIT_FAILURE);
    }
}


// ---------------------------
// Private helper functions
// ---------------------------
void FDGrid::build_sparsity(FullVector<size_type>& ia, FullVector<size_type>& ja) {

    int inode = -1; int icount = -1;

    ia[0] = 0;
    for (size_type iz = 1; iz <= nz; ++iz) {
        for (size_type iy = 1; iy <= ny; ++iy) {
            for (size_type ix = 1; ix <= nx; ++ix) {

                ++inode;
                if (iz > 1) ja[++icount] = inode - nx*ny;
                if (iy > 1) ja[++icount] = inode - nx;
                if (ix > 1) ja[++icount] = inode - 1;

                ja[++icount] = inode;
                if (ix < nx) ja[++icount] = inode + 1;
                if (iy < ny) ja[++icount] = inode + nx;
                if (iz < nz) ja[++icount] = inode + nx*ny;

                ia[inode+1] = icount + 1;
            }
        }
    }

    for (size_type i = 0; i < n; ++i) {
        for (size_type ii = ia[i]; ii < ia[i+1]; ++ii) {
            if (ja[ii == i]) idiag[i] = ii;
        }
    }
}

void FDGrid::set_ascend() {

    FullVector<size_type> ia(n+1);
    FullVector<size_type> ja(7*n);
    FullVector<FloatType>  a(7*n);

    // construct the sparsity pattern
    build_sparsity(ia, ja);

    // fill in the entries of LHS
    for (size_type i = 0; i < n; ++i) {
        int isave = -1; double temp = 0.0;
        for (size_type jj = ia[i]; jj < ia[i+1]; ++jj) {
            if (ja[jj] == i) {
                isave = jj;
            } else {
                a[jj] = -1.0; temp += a[jj];
            }
        }

        if (isave == -1) {
            std::cout << "Error: Row " << i
                << " in ia-ja sparse matrix has a 0 diagonal"
                << std::endl;
            exit(EXIT_FAILURE);
        }
        a[isave] = -temp + 0.1;
    }


    // construct the RHS so that the soln is {1,2,3,4,5...}
    for (size_type i = 0; i < n; ++i) {
        rhs[i] = 0.0;
        for (size_type ii = ia[i]; ii < ia[i+1]; ++ii) {
            rhs[i] += a[ii] * (ja[ii] + 1);
        }
    }

    lhs = SparseMatrixIaja<FloatType>
        (n, std::move(ia), std::move(ja), std::move(a));
    lhs.compress_storage();
}

void FDGrid::set_anisotropic_K() {

    // allocate space
    FullVector<size_type> ia(n+1);
    FullVector<size_type> ja(7*n);
    FullVector<FloatType>  a(7*n);

    // construct the sparsity pattern
    build_sparsity(ia, ja);

    const double bigB = 1.0e8;

    double K, ***Kx, ***Ky, ***Kz;
    size_type i, j, k;
    int isave;
    double temp, val;

    double hx, hy, hz;
    hx = hy = hz = 1.0 / nz;

    // define Kx, Ky, Kz
    Kx = new double** [nx];
    Kx[0] = new double* [nx*ny];
    for (i=1; i < nx; i++) Kx[i] = Kx[i-1]+ny;
    Kx[0][0] = new double[nx*ny*nz];
    for (i=1; i < nx; i++) Kx[i][0] = Kx[i-1][0]+ny*nz;
    for (i=0; i < nx; i++) for (j=1; j < ny; j++) Kx[i][j] = Kx[i][j-1]+nz;

    Ky = new double** [nx];
    Ky[0] = new double* [nx*ny];
    for (i=1; i < nx; i++) Ky[i] = Ky[i-1]+ny;
    Ky[0][0] = new double[nx*ny*nz];
    for (i=1; i < nx; i++) Ky[i][0] = Ky[i-1][0]+ny*nz;
    for (i=0; i < nx; i++) for (j=1; j < ny; j++) Ky[i][j] = Ky[i][j-1]+nz;

    Kz = new double** [nx];
    Kz[0] = new double* [nx*ny];
    for (i=1; i < nx; i++) Kz[i] = Kz[i-1]+ny;
    Kz[0][0] = new double[nx*ny*nz];
    for (i=1; i < nx; i++) Kz[i][0] = Kz[i-1][0]+ny*nz;
    for (i=0; i < nx; i++) for (j=1; j < ny; j++) Kz[i][j] = Kz[i][j-1]+nz;

    for (k=0; k < nz; k++) {
      for (j=0; j < ny; j++) {
        for (i=0; i < nx; i++) {
          Kx[i][j][k] = 1.0;
          if ((i > (nx/4)) && (i < ((3*nx)/4))) {
            if ((j > (ny/4)) && (j < ((3*ny)/4))) {
              if ((k > (nz/4)) && (k < ((3*nz)/4))) {
                Kx[i][j][k] = 1.0;
              }
            }
          }
          Ky[i][j][k] = 1.0;
          if ((i > (nx/4)) && (i < ((3*nx)/4))) {
            if ((j > (ny/4)) && (j < ((3*ny)/4))) {
              if ((k > (nz/4)) && (k < ((3*nz)/4))) {
                Ky[i][j][k] = 1.0;
              }
            }
          }
          Kz[i][j][k] = 10000000.0; // Peter's killer value 10000000.0;
        }
      }
    }

    for (size_type I=0; I < n; I++) {
      isave = -1;
      temp = 0.0;

      k = I/(nx*ny);
      j = (I%(nx*ny))/nx;
      i = (I%(nx*ny))%nx;

      for (size_type ii=ia[I]; ii < ia[I+1]; ii++) {
        size_type col_index = ja[ii];

        // down
        if (col_index == I-nx*ny) {
          K = 2*Kz[i][j][k-1]*Kz[i][j][k]/(Kz[i][j][k-1]+Kz[i][j][k]+1e-10);
          val = -K/(hz*hz);
          a[ii] = val;
          temp -= val;
        }

        // bottom
        if (col_index == I-nx) {
          K = 2*Ky[i][j-1][k]*Ky[i][j][k]/(Ky[i][j-1][k]+Ky[i][j][k]+1e-10);
          val = -K/(hy*hy);
          a[ii] = val;
          temp -= val;
        }

        // left
        if (col_index == I-1) {
          K = 2*Kx[i-1][j][k]*Kx[i][j][k]/(Kx[i-1][j][k]+Kx[i][j][k]+1e-10);
          val = -K/(hx*hx);
          a[ii] = val;
          temp -= val;
        }

        // right
        if (col_index == I+1) {
          K = 2*Kx[i+1][j][k]*Kx[i][j][k]/(Kx[i+1][j][k]+Kx[i][j][k]+1e-10);
          val = -K/(hx*hx);
          a[ii] = val;
          temp -= val;
        }

        // top
        if (col_index == I+nx) {
          K = 2*Ky[i][j+1][k]*Ky[i][j][k]/(Ky[i][j+1][k]+Ky[i][j][k]+1e-10);
          val = -K/(hy*hy);
          a[ii] = val;
          temp -= val;
        }

        // up
        if (col_index == I+nx*ny) {
          K = 2*Kz[i][j][k+1]*Kz[i][j][k]/(Kz[i][j][k+1]+Kz[i][j][k]+1e-10);
          val = -K/(hz*hz);
          a[ii] = val;
          temp -= val;
        }

        // centre
        if (col_index == I) {
          isave = ii;
        }
      }
      assert(isave != -1);// must be diagonal entry
      a[isave] = temp;
    }

    // Kx, Ky, Kz no longer needed
    delete [] Kx[0][0];
    delete [] Kx[0];
    delete [] Kx;
    delete [] Ky[0][0];
    delete [] Ky[0];
    delete [] Ky;
    delete [] Kz[0][0];
    delete [] Kz[0];
    delete [] Kz;

    for (i=0; i < n; i++) {
      rhs[i] = 0.0;
    }
    rhs[0] = -1.0;
    rhs[n-1] = 1.0;

    rhs[n-1] *= bigB;
    a[ idiag[n-1] ] += bigB;

    lhs = SparseMatrixIaja<FloatType>
        (n, std::move(ia), std::move(ja), std::move(a));
    lhs.compress_storage();
}

IAJA_NAMESPACE_CLOSE
