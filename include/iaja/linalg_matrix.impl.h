#include <iaja/iaja_config.h>

#include <cassert>
#include <iomanip>
#include <iostream>
#include <utility>

IAJA_NAMESPACE_OPEN

/* ----------------------
 *  Ctor, Dtor, Assign
 * ---------------------- */

template <typename T>
SparseMatrixIaja<T>::SparseMatrixIaja(size_type nr, size_type nc, size_type nnz,
        const size_type* ia, const size_type* ja):
    SparsityIaja(nr, nc, nnz, ia, ja), a(nnz) {}


template <typename T>
SparseMatrixIaja<T>::SparseMatrixIaja(size_type nr, size_type nc, size_type nnz,
        const size_type* ia, const size_type* ja, const T& a):
    SparsityIaja(nr, nc, nnz, ia, ja), a(nnz, a) {}


template <typename T>
SparseMatrixIaja<T>::SparseMatrixIaja(size_type nr, size_type nc,  size_type nnz,
        const size_type* ia, const size_type* ja, const T* a):
    SparsityIaja(nr, nc, nnz, ia, ja), a(nnz, a) {}

template <typename T>
SparseMatrixIaja<T>::SparseMatrixIaja(size_type nc,
        const FullVector<size_type>& ia,
        const FullVector<size_type>& ja,
        const FullVector<T>& a):
    SparsityIaja(nc, ia, ja), a(a) 
{ assert(ja.length() == a.length()); }


template <typename T>
SparseMatrixIaja<T>::SparseMatrixIaja(size_type nc,
        FullVector<size_type>&& ia,
        FullVector<size_type>&& ja,
        FullVector<T>&& a):
    SparsityIaja(nc, std::move(ia), std::move(ja)),
    a(std::move(a))
{ assert(ja.length() == a.length()); }


template <typename T>
SparseMatrixIaja<T>::SparseMatrixIaja(size_type nc,
        const std::vector<size_type>& ia,
        const std::vector<size_type>& ja,
        const std::vector<T>& a):
    SparsityIaja(nc, ia, ja), a(a)
{ assert(ja.size() == a.size()); }


template <typename T>
SparseMatrixIaja<T>::SparseMatrixIaja(SparseMatrixIaja<T>&& rhs):
    SparsityIaja(std::move(rhs)), a(std::move(rhs.a)) {}


template <typename T>
SparseMatrixIaja<T>& SparseMatrixIaja<T>:: operator=(SparseMatrixIaja<T>&& rhs) {
    if (this != &rhs) {
        SparsityIaja::operator=(std::move(rhs));
        a = std::move(rhs.a);
    }
    return *this;
}

/* ----------
 *    I/O
 * ---------- */

template <typename T>
std::ostream& operator<<(std::ostream& os, const SparseMatrixIaja<T>& mtrx) {

    using size_type = typename SparseMatrixIaja<T>::size_type;

    for (size_type i = 0, ctr = 0; i < mtrx.nr; ++i) {
        for (size_type j = 0; j < mtrx.nc ; ++j) {
            if ( ctr < mtrx.ia[i+1] && mtrx.ja[ctr] == j ) {
                os << mtrx.a[ctr++] << " ";
            } else {
                os << std::right << 0 << " ";
            }
        } // column loop
        os << "\n";
    }  // row loop

    return os;
}


template <> inline
std::ostream& operator<<(std::ostream& os, const SparseMatrixIaja<double>& mtrx) {

    using size_type = SparseMatrixIaja<double>::size_type;
    auto width = PRINT_WIDTH_DOUBLE;

    for (size_type i = 0, ctr = 0; i < mtrx.nr; ++i) {
        for (size_type j = 0; j < mtrx.nc ; ++j) {
            if ( ctr < mtrx.ia[i+1] && mtrx.ja[ctr] == j ) {
                os << std::setw(width) << std::scientific
                    << std::setprecision(PRINT_PRECISION_DOUBLE)
                    << std::right << mtrx.a[ctr++];
            } else {
                os << std::setw(width) << std::scientific << std::right << 0;
            }
        } // column loop
        os << "\n";
    }  // row loop

    return os;
}

template <typename T>
std::ostream& SparseMatrixIaja<T>::print_compressed(std::ostream& os) const {
    os << nr <<"\n" << nnz << "\n" <<
        ia << "\n" << ja << "\n" << a << "\n";
    return os;
}

/* ---------------
 *    INDEXING
 * --------------- */

template <typename T>
T& SparseMatrixIaja<T>:: operator[](size_type i) {
    assert(i >= 0 && i < nnz);
    return a[i];
}

template <typename T>
const T& SparseMatrixIaja<T>:: operator[](size_type i) const {
    assert(i >= 0 && i < nnz);
    return a[i];
}


/* ----------------
 *    NUMERICS
 * ---------------- */

template <typename T>
FullVector<T> SparseMatrixIaja<T>:: operator* (const FullVector<T>& x) const {
    assert(nc == x.length());
    FullVector<T> y(nr, T());
    for (size_type i = 0; i < nr; ++i) {
        for (size_type jj = ia[i]; jj < ia[i+1]; ++jj) {
            y[i] += a[jj] * x[ja[jj]];
        }
    }
    return y;
}


template <typename T>
SparseMatrixIaja<T> SparseMatrixIaja<T>::transpose() const {

    // allocate transposed matrix
    SparseMatrixIaja<T> At(nc, nr, nnz);

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

template <typename T>
void SparseMatrixIaja<T>::compress_storage() {
    SparsityIaja::compress_storage();
    a.resize(nnz);
}

IAJA_NAMESPACE_CLOSE
