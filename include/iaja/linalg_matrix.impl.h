#include <iaja/global_defs.h>

#include <cassert>
#include <iomanip>

IAJA_NAMESPACE_OPEN

/* -------------------------------------------
 * Constructors and destructors 
 * ------------------------------------------- */
template <typename T>
SparseMatrix<T>::SparseMatrix(unsigned int nr, unsigned int nc,  unsigned int nnz,
        const unsigned int* ia, const unsigned int* ja, const double* a):
    ia(nr+1, ia), ja(nnz, ja), a(nnz, a), nr(nr), nc(nc), nnz(nnz) {}

template <typename T>
SparseMatrix<T>::SparseMatrix(unsigned int nr, unsigned int nc, unsigned int nnz)
    : ia(nr+1), ja(nnz), a(nnz), nr(nr), nc(nc), nnz(nnz) {}


template <typename T>
SparseMatrix<T>::SparseMatrix(const SparseMatrix<T>& mat)
    : nr(mat.nr), nc(mat.nc), nnz(mat.nnz), ia(mat.ia), ja(mat.ja), a(mat.a) {}


/* -------------------------------------------
 * Operator overloading
 * ------------------------------------------- */
template <typename T>
SparseMatrix<T>& SparseMatrix<T>:: operator= (const SparseMatrix<T>& rhs) {
    if (this != &rhs) {
        ja = rhs.ja; ia = rhs.ia; a = rhs.a;
        nr = rhs.nr; nc = rhs.nc; nnz = rhs.nnz;
    }
    return *this;
}

template <typename T>
std::ostream& operator<< (std::ostream& os, const SparseMatrix<T>& mtrx) {

    const unsigned int width = PRINT_WIDTH_DOUBLE;
    unsigned int ctr = 0;

    // row loop
    for (unsigned int i = 0; i < mtrx.nr; ++i) {

        // column loop
        for (unsigned int j = 0; j < mtrx.nc ; ++j) {
            if (mtrx.ja[ctr] == j) {
                os << std::setw(width) << std::scientific << std::right << mtrx.a[ctr] << " ";
                ++ctr;
            } else {
                os << std::setw(width) << std::scientific << std::right << 0 << " ";
            }
        } // j
        os << "\n";
    }  // i

    return os;
}

template <typename T>
T& SparseMatrix<T>:: operator[](unsigned int i) {
    assert(i >= 0 && i < nnz);
    return a[i];
}

template <typename T>
const T& SparseMatrix<T>:: operator[](unsigned int i) const {
    assert(i >= 0 && i < nnz);
    return a[i];
}

template <typename T>
FullVector<T> SparseMatrix<T>:: operator* (const FullVector<T> x) const {

    assert(nc == x.length());
    FullVector<T> y(nr);

    for (unsigned int i = 0; i < nc; ++i) {
        for (unsigned int jj = ia[i]; jj < ia[i+1]; ++jj) {
            y[i] += a[jj] * x[ja[jj]];
        }
    }
    return y;
}

IAJA_NAMESPACE_CLOSE
