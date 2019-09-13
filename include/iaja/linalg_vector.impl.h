#include <iaja/iaja_config.h>
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

/* ------ Constructors and Destructors ------ */
template <typename T>
FullVector<T>::FullVector() {
    n = 0; a = nullptr;
}

template <typename T>
FullVector<T>::FullVector(size_type n_in) {
    // allocate and 0-init
    n = n_in;
    a = new T[n]();
}

template <typename T>
FullVector<T>::FullVector(size_type n_in, const T* a_in) {
    n = n_in;
    a = new T[n_in];
    for (size_type i = 0; i < n; ++i)
        a[i] = a_in[i];
}

template <typename T>
FullVector<T>::FullVector(const FullVector<T>& v) {
    // copy constructor will allocate new space
    n = v.n;
    a = new T[n];
    for (size_type i = 0; i < n; ++i)
        a[i] = v[i];
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


/* ------ Operator Overloading ------ */
template <typename T>
const T& FullVector<T>:: operator[](size_type i) const {
    assert(i >= 0 && i < n);
    return a[i];
    // This also works... but tedious
    // return (*(const_cast<FullVector<T>*>(this)))[i];
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
            a[ctr] = *i; ++ctr;
        }
    }
    return *this;
}

template <typename T>
std::ostream& operator<< (std::ostream& os, const FullVector<T>& x) {
    const unsigned int width = PRINT_WIDTH_DOUBLE;
    for (decltype(x.n) i = 0; i < x.n; ++i)
            os << std::setw(width) << std::scientific <<std::right<< x.a[i];
    return os;
}


/* ------ Public Methods ------- */
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
void FullVector<T>::multiply(const SparseMatrix<T>& A, const FullVector<T>& x) {

    assert(A.ncol() == x.length());
    assert(A.nrow() == this->length());

    for (size_type i = 0; i < A.nrow(); ++i) {
        a[i] = 0.0;
        for (size_type jj = A.ia[i]; jj < A.ia[i+1]; ++jj) {
            a[i] += A[jj] * x[ A.ja[jj] ];
        }
    }
}

template <typename T>
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

/* ------ Constructors and Destructors ------ */
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
}

template <typename T>
SparseVector<T>::SparseVector(const SparseVector& rhs)
    :n(rhs.n), nnz(rhs.nnz), ja(rhs.ja), a(rhs.a) {
}

template <typename T>
SparseVector<T>::SparseVector(SparseVector&& rhs) :
    n(rhs.n), nnz(rhs.nnz), ja(std::move(rhs.ja)), a(std::move(rhs.a)) {
    rhs.n = rhs.nnz = 0;
}


/* ------ Operator Overloading ------ */
template <typename T>
std::ostream& operator<<(std::ostream& os, const SparseVector<T>& vin) {

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

IAJA_NAMESPACE_CLOSE
