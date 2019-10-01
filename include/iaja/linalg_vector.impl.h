
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

/* ---------------------
 *   Ctor, Dtor, Assign 
 * --------------------- */

template <typename T> inline
FullVector<T>::FullVector(size_type n, const T& val):
    n(n), a(new T[n]) {
    for (size_type i = 0; i < n; ++i ) a[i] = val;
}

template <typename T> inline
FullVector<T>::FullVector(size_type n, const T* arr):
    n(n), a(new T[n]) {
    for (size_type i = 0; i < n; ++i) a[i] = arr[i];

#ifdef SHOW_CALLED_FCN
    std::cout << "FullVector pointer constructor called\n";
#endif
}

template <typename T> inline
FullVector<T>::FullVector(size_type n, T* && arr):
    n(n), a(arr) {
    arr = nullptr;

#ifdef SHOW_MOVE_COPY
    std::cout << "FullVector rvalue ref constructor called\n";
#endif
}

template <typename T> inline
FullVector<T>::FullVector(const std::vector<T>& rhs):
    n(rhs.size()), a(new T[rhs.size()]) {
    size_type ctr = 0;
    for (auto i = rhs.cbegin(); i != rhs.cend(); ++i) {
        a[ctr++] = *i;
    }
#ifdef SHOW_MOVE_COPY
    std::cout << "FullVector std::vector constructor called\n";
#endif
}

template <typename T>
FullVector<T>::FullVector(const FullVector<T>& v)
    : n(v.n), a(new T[n]) {
    for (size_type i = 0; i < n; ++i) a[i] = v[i];
#ifdef SHOW_MOVE_COPY
    std::cout << "FullVector copy constructor called\n";
#endif
}

template <typename T>
FullVector<T>::FullVector(FullVector<T>&& v)
    : n(v.n), a(v.a) {
    v.n = 0; v.a = nullptr;
#ifdef SHOW_MOVE_COPY
    std::cout << "FullVector move constructor called\n";
#endif
}

template <typename T>
FullVector<T>::~FullVector() {
    delete[] a; a = nullptr;
}

template <typename T>
FullVector<T>& FullVector<T>:: operator= (const FullVector<T>& rhs) {
    if (this != &rhs) {
        if (n != rhs.n) {
            delete[] a; // delete nullptr is safe
            a = new T[rhs.n];
        }
        n = rhs.n;
        for (size_type i = 0; i < n; ++i) a[i] = rhs[i];
    }

#ifdef SHOW_MOVE_COPY
    std::cout << "FullVector copy assign called\n";
#endif
    return *this;
}

template <typename T>
FullVector<T>& FullVector<T>:: operator= (FullVector<T>&& rhs) {
    if (this != &rhs) {
        delete[] a;
        a = rhs.a; n = rhs.n;
        rhs.a = nullptr; rhs.n = 0;
    }

#ifdef SHOW_MOVE_COPY
    std::cout << "FullVector move assign called\n";
#endif
    return *this;
}

/* ------------
 *  Indexing
 * ------------ */

template <typename T> inline
const T& FullVector<T>:: operator[](size_type i) const {
    assert(i >= 0 && i < n);
    // return (*(const_cast<FullVector<T>*>(this)))[i];
    return a[i];
}

template <typename T> inline
T& FullVector<T>:: operator[](size_type i) {
    assert(i >= 0 && i < n);
    return a[i];
}

template <typename T> inline
T& FullVector<T>:: operator[](long int i) {
    size_type ii = (i < 0) ? ( n + i ) : i;
    assert(ii >= 0 && static_cast<size_type>(i) < n);
    return a[ii];
}

template <typename T> inline
const T& FullVector<T>::operator[](long int i) const {
    // reusing [] defined for non-const
    return (*const_cast< FullVector<T>* >(this))[i];
}

template <typename T> inline
T& FullVector<T>:: operator[](int i) {
    int ii = (i < 0) ? ( n + i ) : i;
    assert(ii >= 0 && static_cast<size_type>(ii) < n);
    return a[ii];
}

template <typename T> inline
const T& FullVector<T>:: operator[](int i) const {
    return ( *const_cast< FullVector<T>* > (*this) )[i];
}


/* -------------------
 *  Size / Dimension
 * ------------------- */

template <typename T>
void FullVector<T>::resize(size_t new_leng) {
    if (new_leng != n) {
        if (new_leng == 0) {
            delete[] a; a = nullptr; n = 0;
        }
        T* a_resz = new T[new_leng]();
        n = (new_leng > n) ? n : new_leng;
        for (size_type i = 0; i < n; ++i) {a_resz[i] = a[i];}
        delete[] a;
        a = a_resz; a_resz = nullptr;
    }
}

template <typename T>
void FullVector<T>::remove_trailing_zeros() {
    unsigned int num2delete = 0;
    auto i = cend();
    while ( *(--i) == T() && i != cbegin() ) { ++num2delete; }
    resize(n - num2delete);
}

/* ---------
 *   I/O 
 * --------- */

template <typename T> inline
std::ostream& operator<< (std::ostream& os, const FullVector<T>& x) {
    for (auto i = x.cbegin(); i != x.cend(); ++i)
        os << *i << " ";
    return os;
}

template <> inline
std::ostream& operator<< (std::ostream& os, const FullVector<double>& x) {
    const unsigned int width = PRINT_WIDTH_DOUBLE;
    for (auto i = x.cbegin(); i != x.cend(); ++i)
            os << std::setw(width) << std::scientific
                << std::setprecision(PRINT_PRECISION_DOUBLE)
                << std::left << *i;
    return os;
}


/* -------------------------
 *   Numerical Operations 
 * ------------------------- */

// ------------------------------------------
// the following routines will put the result
// into the current calling object
// ------------------------------------------
//
template <typename T> inline
void FullVector<T>::cumsum() {
    for ( size_type i = 1; i < n; ++i )
        a[i] += a[i-1];
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
        a[i] = T();
        for (size_type jj = A.ia[i]; jj < A.ia[i+1]; ++jj) {
            a[i] += A[jj] * x[ A.ja[jj] ];
        }
    }
}

template <typename T> inline
void FullVector<T>::multiply(const T& s) {
    for ( size_type i = 0; i < n; ++i ) a[i] = s*a[i];
}

template <typename T> inline
void FullVector<T>::saxpy(T alpha,
        const FullVector<T>& x,
        const FullVector<T>& y) {
    // a = alpha*x + y
    for (size_type i = 0; i < n; i++)
        a[i] = alpha*x[i] + y[i];
}


// ---------------------------------------
// these methods will return the result
// as a newly constructed object
// ---------------------------------------

template <typename T> inline
T FullVector<T>:: operator *(const FullVector<T>& x) const {
    assert(n == x.n);
    T s = T();
    for (size_type i = 0 ; i < x.n; ++i)
        s += a[i] * x[i];
    return s;
}

template <typename T> inline
FullVector<T>& FullVector<T>:: operator *=(const T& s) {
    for (auto i = begin(); i != end(); ++i) *i = *i * s;
    return *this;
}

template <typename T> inline
FullVector<T>&& FullVector<T>:: operator +(FullVector<T>&& other) const {
    assert(n == other.length());
    for (size_t i = 0; i < n; ++i)
        other[i] = a[i] + other[i];

#ifdef SHOW_CALLED_FCN
    std::cout << "operator+ with rvalue ref called\n";
#endif

    return std::move(other);
}

template <typename T> inline
FullVector<T> FullVector<T>:: operator +(const FullVector<T>& other) const {
    assert(n == other.length());
    FullVector<T> result(n);
    for (size_t i = 0; i < n; ++i)
        result[i] = a[i] + other[i];

#ifdef SHOW_CALLED_FCN
    std::cout << "operator+ with lvalue ref called\n";
#endif

    return result;
}

template <typename T> inline
FullVector<T> FullVector<T>:: operator -(const FullVector<T>& other) const {
    assert(n == other.length());
    FullVector<T> result(n);
    for (size_t i = 0; i < n; ++i)
        result[i] = a[i] - other[i];
    return result;
}

template <typename T> inline
T FullVector<T>::norm_l2() const {
    return sqrt( (*this) * (*this) );
}

// TODO infinity norm

/* ============================================ *
 *                SPARSEVECTOR                  *
 * ============================================ */

/* ---------------------
 *  Ctor, Dtor, Assign 
 * --------------------- */

template <typename T>
SparseVector<T>::SparseVector():
    ja(), a() {
    n = nnz = 0;
}

template <typename T>
SparseVector<T>::SparseVector(size_type n, size_type nnz, const size_type* ja):
    n(n), nnz(nnz), ja(nnz, ja), a(nnz) {}

template <typename T>
SparseVector<T>::SparseVector(size_type n, size_type nnz,
        const size_type* ja, const T& val):
    n(n), nnz(nnz), ja(nnz, ja), a(nnz, val) {}

template <typename T>
SparseVector<T>::SparseVector(size_type n_in, size_type nnz_in,
        const size_type* ja_in, const T* a_in):
    n(n_in), nnz(nnz_in), ja(nnz_in, ja_in), a(nnz_in, a_in) {
    assert(nnz_in <= n_in);
    for (size_type i = 0; i < nnz_in-1; ++i)
        assert(ja_in[i] < ja_in[i+1]);
}

template <typename T>
SparseVector<T>::SparseVector(size_type n,
        const FullVector<size_type>& ja, const FullVector<T>& a):
    n(n), nnz(a.length()), ja(ja), a(a)
{ assert(ja.length() == a.length()); }


template <typename T>
SparseVector<T>::SparseVector(size_type n,
        FullVector<size_type>&& ja, FullVector<T>&& a):
    n(n), nnz(a.length()),
    ja(std::move(ja)), a(std::move(a))
{ assert(ja.length() == a.length()); }

template <typename T>
SparseVector<T>::SparseVector(size_type n, FullVector<size_type>&& ja):
    n(n), nnz(ja.length()), ja(std::move(ja)), a(nnz) {}


template <typename T>
SparseVector<T>::SparseVector(size_type n,
        const std::vector<size_type>& ja,
        const std::vector<T>& a) :
    n(n), nnz(a.size()), ja(ja), a(a)
{ assert(ja.size() == a.size()); }


template <typename T>
SparseVector<T>::SparseVector(SparseVector&& rhs)
    : n(rhs.n), nnz(rhs.nnz),
    ja(std::move(rhs.ja)), a(std::move(rhs.a)) {
    rhs.n = rhs.nnz = 0;
}


template <typename T>
SparseVector<T>& SparseVector<T>:: operator=(SparseVector<T>&& rhs) {
    if (this != &rhs) {
        n = rhs.n; nnz = rhs.nnz;
        rhs.nnz = rhs.n = 0;
        ja = std::move(rhs.ja);
        a = std::move(rhs.a);
    }
    return *this;
}

/* ---------------
 *    Indexing 
 * --------------- */

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

/* ------
 *  I/O
 * ------ */

template <typename T>
std::ostream& operator<<(std::ostream& os, const SparseVector<T>& vin) {
    for (decltype(vin.n) i = 0, j = 0; i < vin.n; ++i) {
        if (j < vin.nnz && i == vin.ja[j]) {
            os << vin.a[j++] << " ";
        } else {
            os << 0 << " ";
        }
    }
    return os;
}

template <> inline
std::ostream& operator<<(std::ostream& os, const SparseVector<double>& vin) {

    const unsigned int width = PRINT_WIDTH_DOUBLE;

    for (decltype(vin.n) i = 0, j = 0; i < vin.n; ++i) {
        if (j < vin.nnz && i == vin.ja[j]) {
            os << std::setw(width) << std::scientific <<std::left<< vin.a[j++];
        } else {
            os << std::setw(width) << std::scientific << std::left << 0;
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
