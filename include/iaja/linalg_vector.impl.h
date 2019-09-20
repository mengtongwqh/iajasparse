
#include <iaja/linalg_vector.h>

#include <cassert>
#include <cmath>
#include <iomanip>
#include <vector>
#include <iostream>
#include <utility>

IAJA_NAMESPACE_OPEN

/* ============================================ *
 *                  VECTOR                      *
 * ============================================ */

/* ------ Ctor, Dtor, Assign ------ */
template <typename T>
FullVector<T>::FullVector(size_type n, const T* arr)
    : n(n), a(new T[n]) {
    for (size_type i = 0; i < n; ++i)
        a[i] = arr[i];
}

template <typename T>
FullVector<T>::FullVector(const FullVector<T>& v)
    : n(v.n), a(new T[n]) {
    for (size_type i = 0; i < n; ++i) a[i] = v[i];
}

template <typename T>
FullVector<T>::FullVector(FullVector<T>&& v)
    : n(v.n), a(v.a) {
    v.n = 0; v.a = nullptr;
}

template <typename T>
FullVector<T>::~FullVector() {
    delete[] a; a = nullptr;
}

template <typename T>
FullVector<T>& FullVector<T>:: operator= (const FullVector<T>& rhs) {
    if (this != &rhs) {
        if (n != rhs.n) {
            if (n != 0) delete[] a;
            a = new T[rhs.n];
        }
        n = rhs.n;
        for (size_type i = 0; i < n; ++i) a[i] = rhs[i];
    }
    return *this;
}

template <typename T>
FullVector<T>& FullVector<T>:: operator= (const std::vector<T>& rhs) {
    if (n != rhs.size()) {
        if (n != 0) {
            delete[] a;
            a = new T[rhs.size()];
        }
        n = rhs.size();
        size_type ctr = 0;
        for (auto i = rhs.cbegin(); i != rhs.cend(); ++i) {
            a[ctr++] = *i;
        }
    }
    return *this;
}

template <typename T>
FullVector<T>& FullVector<T>:: operator= (FullVector<T>&& rhs) {
    if (this != &rhs) {
        if (n != 0) delete[] a;
        a = rhs.a; n = rhs.n;
        rhs.a = nullptr; rhs.n = 0;
    }    
    return *this;
}

/* ------ Indexing ------ */
template <typename T>
const T& FullVector<T>:: operator[](size_type i) const {
    assert(i >= 0 && i < n);
    // return (*(const_cast<FullVector<T>*>(this)))[i];
    return a[i];
}

template <typename T>
T& FullVector<T>:: operator[](size_type i) {
    assert(i >= 0 && i < n);
    return a[i];
}

template <typename T>
T& FullVector<T>:: operator[](long int i) {
    size_type ii = (i < 0) ? ( n + i ) : i;
    assert(ii >= 0 && static_cast<size_type>(i) < n);
    return a[ii];
}

template <typename T>
const T& FullVector<T>::operator[](long int i) const {
    // reusing [] defined for non-const
    return (*const_cast< FullVector<T>* >(this))[i];
}

template <typename T>
T& FullVector<T>:: operator[](int i) {
    int ii = (i < 0) ? ( n + i ) : i;
    assert(ii >= 0 && static_cast<size_type>(ii) < n);
    return a[ii];
}

template <typename T>
const T& FullVector<T>:: operator[](int i) const {
    return ( *const_cast< FullVector<T>* > (*this) )[i];
}

template <typename T>
T FullVector<T>:: operator* (const FullVector<T>& x) const {
    assert(n == x.n);
    T s = 0.0;
    for (size_type i = 0 ; i < x.n; ++i)
        s += a[i] * x[i];
    return s;
}


/* ----- I/O ----- */
template <typename T>
std::ostream& operator<< (std::ostream& os, const FullVector<T>& x) {
    for (decltype(x.n) i = 0; i < x.n; ++i) os << x.a[i] << " "; 
    return os;
}

template <> inline
std::ostream& operator<< (std::ostream& os, const FullVector<double>& x) {
    const unsigned int width = PRINT_WIDTH_DOUBLE;
    for (decltype(x.n) i = 0; i < x.n; ++i)
            os << std::setw(width) << std::scientific <<std::right<< x.a[i];
    return os;
}


/* ------ Numerical Operations ------- */
template <typename T> inline
FullVector<T>& FullVector<T>:: operator*=(const T& s) {
    for (auto i = begin(); i != end(); ++i) *i = *i * s;
    return *this;
}
 
template <typename T>
T FullVector<T>::norm_l2() const {
    return sqrt( (*this) * (*this) );
}

template <typename T>
void FullVector<T>::saxpy(T alpha, const FullVector<T>& x, const FullVector<T>& y) {
    // a = alpha*x + y
    for (size_type i = 0; i < n; i++)
        a[i] = alpha*x[i] + y[i];
}

template <typename T>
void FullVector<T>::add(const FullVector<T>& b, const FullVector<T>& c) {
    // a = b - c
    assert(b.length() == c.length());
    for (size_type i = 0; i < n; i++)
        a[i] = b[i] + c[i];
}

template <typename T>
void FullVector<T>::minus(const FullVector<T>& b, const FullVector<T>& c) {
    // a = b - c
    assert(b.length() == c.length());
    for (size_type i = 0; i < n; i++)
        a[i] = b[i] - c[i];
}

template <typename T>
void FullVector<T>::multiply(const SparseMatrixIaja<T>& A, const FullVector<T>& x) {

    assert(A.ncol() == x.length());
    assert(A.nrow() == this->length());

    for (size_type i = 0; i < A.nrow(); ++i) {
        a[i] = 0.0;
        for (size_type jj = A.ia[i]; jj < A.ia[i+1]; ++jj) {
            a[i] += A[jj] * x[ A.ja[jj] ];
        }
    }
}

template <typename T> inline
void FullVector<T>::multiply(const T& s) {
    for ( size_type i = 0; i < n; ++i ) a[i] = s*a[i];
}


template <typename T>
void FullVector<T>::cumsum() {
    for ( size_type i = 1; i < n; ++i )
        a[i] += a[i-1];
}

/* ============================================ *
 *                SPARSEVECTOR                  *
 * ============================================ */

/* ------ Ctor, Dtor, Assign ------ */
template <typename T>
SparseVector<T>::SparseVector():
    ja(), a() {
    n = nnz = 0;
}

template <typename T>
SparseVector<T>::SparseVector(size_type n_in, size_type nnz_in):
    n(n_in), nnz(nnz_in), ja(nnz_in), a(nnz_in) {
    assert(nnz_in <= n_in);
}

template <typename T>
SparseVector<T>::SparseVector(size_type n_in, size_type nnz_in,
        const size_type* ja_in, const T* a_in):
    n(n_in), nnz(nnz_in), ja(nnz_in, ja_in), a(nnz_in, a_in) {
    assert(nnz_in <= n_in);
    for (size_type i = 0; i < nnz_in-1; ++i)
        assert(ja_in[i] < ja_in[i+1]);
}

template <typename T>
SparseVector<T>::SparseVector(const SparseVector& rhs)
    :n(rhs.n), nnz(rhs.nnz), ja(rhs.ja), a(rhs.a) {
}

template <typename T>
SparseVector<T>::SparseVector(SparseVector&& rhs)
    : n(rhs.n), nnz(rhs.nnz),
    ja(std::move(rhs.ja)), a(std::move(rhs.a)) {
    rhs.n = rhs.nnz = 0;
}

template <typename T>
SparseVector<T>& SparseVector<T>:: operator=(const SparseVector<T>& rhs) {
    if (this != &rhs) {
        n = rhs.n; nnz = rhs.nnz;
        ja = rhs.ja; a = rhs.a;
    }
    return *this;
}

template <typename T>
SparseVector<T>& SparseVector<T>:: operator=(SparseVector<T>&& rhs) {
    if (this != rhs) {
        n = rhs.n; nnz = rhs.nnz;
        rhs.nnz = rhs.n = 0;
        ja = std::move(rhs.ja);
        a = std::move(rhs.a);
    }
    return *this;
}

/* ----- Indexing ----- */
template <typename T>
T& SparseVector<T>:: operator[](size_type i) {
    assert(i >= 0 && i < nnz);
    return a[i];
}

template <typename T>
const T& SparseVector<T>:: operator[](size_type i) const {
    assert(i >= 0 && i < nnz);
    return a[i];
}

/* ------ I/O ------ */
template <typename T>
std::ostream& operator<<(std::ostream& os, const SparseVector<T>& vin) {

    for (decltype(vin.n) i = 0, j = 0; i < vin.n; ++i) {
        if (j < vin.nnz && i == vin.ja[j]) {
            os << vin.a[j] << " ";
            ++j;
        } else {
            os << 0 << " ";
    }

    return os;
}

template <> inline
std::ostream& operator<<(std::ostream& os, const SparseVector<double>& vin) {

    const unsigned int width = PRINT_WIDTH_DOUBLE;

    for (decltype(vin.n) i = 0, j = 0; i < vin.n; ++i) {
        if (j < vin.nnz && i == vin.ja[j]) {
            os << std::setw(width) << std::scientific <<std::right<< vin.a[j];
            ++j;
        } else {
            os << std::setw(width) << std::scientific << std::right << 0;
        }
    }

    return os;
}

template <typename T> inline
std::ostream& SparseVector<T>:: print_compressed(std::ostream& os) const {
    os << ja << "\n" << a << "\n";
    return os;
}


IAJA_NAMESPACE_CLOSE
