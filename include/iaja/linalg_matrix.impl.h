#include <iaja/iaja_config.h>

#include <cassert>
#include <iomanip>

IAJA_NAMESPACE_OPEN

/* -------------------------------------------
 * Constructors and destructors 
 * ------------------------------------------- */
template <typename T>
SparseMatrix<T>::SparseMatrix(size_type nr, size_type nc,  size_type nnz,
        const size_type* ia, const size_type* ja, const double* a):
    ia(nr+1, ia), ja(nnz, ja), a(nnz, a), nr(nr), nc(nc), nnz(nnz) {}

template <typename T>
SparseMatrix<T>::SparseMatrix(size_type nr, size_type nc, size_type nnz)
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

    using size_type = decltype(mtrx.nr);
    const unsigned int width = PRINT_WIDTH_DOUBLE;
    size_type ctr = 0;

    // row loop
    for (size_type i = 0; i < mtrx.nr; ++i) {

        // column loop
        for (size_type j = 0; j < mtrx.nc ; ++j) {

            if ( ctr < mtrx.ja.length() && mtrx.ja[ctr] == j ) {
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
T& SparseMatrix<T>:: operator[](size_type i) {
    assert(i >= 0 && i < nnz);
    return a[i];
}

template <typename T>
const T& SparseMatrix<T>:: operator[](size_type i) const {
    assert(i >= 0 && i < nnz);
    return a[i];
}

template <typename T>
FullVector<T> SparseMatrix<T>:: operator* (const FullVector<T> x) const {

    assert(nc == x.length());
    FullVector<T> y(nr);

    for (size_type i = 0; i < nc; ++i) {
        for (size_type jj = ia[i]; jj < ia[i+1]; ++jj) {
            y[i] += a[jj] * x[ja[jj]];
        }
    }
    return y;
}

/* -------------------------------------------
 * Methods
 * ------------------------------------------- */
template <typename T>
SparseMatrix<T> SparseMatrix<T>::transpose() const {

    // allocate transposed matrix
    SparseMatrix<T> At(nc, nr, nnz);

    // row indices
    for (size_type i = 0; i < nnz; ++i)
        ++( At.ia[ja[i]+1] );
    At.ia.cumsum();

    assert(At.ia[nc] == nnz);
    assert(At.ia[0] == 0);

    // column position for each row
    FullVector<size_type> col_idx(At.ia);

    for (size_type i = 0; i < nr; ++i) {
        for (size_type jj = ia[i]; jj < ia[i+1]; ++jj) {
            size_type p = col_idx[ja[jj]]++;
            At.ja[p] = i;
            At.a[p]  = a[jj];
        }
    }

    return At;
}

IAJA_NAMESPACE_CLOSE
