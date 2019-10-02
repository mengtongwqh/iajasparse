#ifndef _LINALG_VECTOR_H_
#define _LINALG_VECTOR_H_

#include <iaja/iaja_config.h>

#include <iostream>
#include <vector>

IAJA_NAMESPACE_OPEN

template <typename T> class SparseMatrixIaja;

/* A little more than the a naked pointer to a block of heap */
template<typename T>
class FullVector {

 public:
  /* Type Aliases */
  using value_type      = T;
  using size_type       = SizeType;

  using pointer         = value_type *;
  using const_pointer   = const value_type *;

  using iterator        = value_type *;
  using const_iterator  = const value_type *;

  using reference       = value_type &;
  using const_reference = const value_type &;

  /* ctor */
  FullVector() : n(0), a(nullptr) {}

  explicit FullVector(size_type n) : n(n), a(new T[n]) {}

  FullVector(size_type n, const T& a);

  FullVector(size_type n, const T* a);

  explicit FullVector(const std::vector<T>& rhs);

  FullVector(size_type n, T* && a);

  /* copy ctor */
  FullVector(const FullVector<T>& v);
  /* copy assign */
  FullVector<T>& operator= (const FullVector<T>& rhs);

  /* move ctor */
  FullVector(FullVector<T>&& v);
  /* move assign */
  FullVector<T>& operator= (FullVector<T>&& rhs);
  /* dtor */
  virtual ~FullVector();

  /* Iterator */
  iterator begin() { return a; }
  iterator end() { return a+n; }

  const_iterator cbegin() const { return a; }
  const_iterator cend() const { return a+n; }

  /* Conversion Operator */
  explicit operator T*() { return a; } // collapse to pointer to array

  /* Indexing */
  reference operator[] (size_type i);
  const_reference operator[] (size_type i) const;

  reference operator[] (long int i);
  const_reference operator[] (long int i) const;

  reference operator[] (int i);
  const_reference operator[] (int i) const;

  /* I/O */
  template <typename U>
  friend std::ostream& operator<< (std::ostream& os, const FullVector<U>& x);

  /* Size & Dimension */
  size_type length() const {return n;}

  void resize(size_t new_leng);

  void remove_trailing_zeros();

  /* Numerical Operation, results returned as objects */
  value_type norm_l2() const;
  value_type operator * (const FullVector<T>& x) const; // dot product

  FullVector<T>   operator + (const FullVector<T>& other) const;
  FullVector<T>&& operator + (FullVector<T>&& other) const;

  FullVector<T>  operator - (const FullVector<T>& other) const; // vector add
  FullVector<T>& operator*= (const T& s); // scalar multiply

  /* Numerical Operation, inplace */
  void cumsum();

  void add(const FullVector<T>& b, const FullVector<T>& c);

  void minus(const FullVector<T>& b, const FullVector<T>& c);

  void multiply(const SparseMatrixIaja<T>& A, const FullVector<T>& x);
  void multiply(const T& s);

  void saxpy(T alpha, const FullVector<T>& x, const FullVector<T>& y);

 protected:
  size_type n;
  T* a;
};


template <typename T>
class SparseVector {

 public:
  using value_type      = T;
  using size_type       = SizeType;

  using pointer         = value_type *;
  using const_pointer   = const value_type *;

  using iterator        = value_type *;
  using const_iterator  = const value_type *;

  using reference       = value_type &;
  using const_reference = const value_type &;

  /* ctor */
  SparseVector();

  SparseVector(size_type n, size_type nnz, const size_type* ja);

  SparseVector(size_type n, size_type nnz, const size_type* ja, const T& val);

  SparseVector(size_type n, size_type nnz, const size_type* ja, const T* a);

  SparseVector(size_type n, const FullVector<size_type>& ja, const FullVector<T>& a);

  SparseVector(size_type n, FullVector<size_type>&& ja);

  SparseVector(size_type n, FullVector<size_type>&& ja, FullVector<T>&& a);

  SparseVector(size_type n, const std::vector<size_type>& ja, const std::vector<T>& a);

  /* copy ctor */
  SparseVector(const SparseVector<T>& rhs) = default;
  /* copy assign */
  SparseVector<T>& operator= (const SparseVector<T>& rhs) = default;

  /* move ctor */
  SparseVector(SparseVector<T>&& rhs);
  /* move assign */
  SparseVector<T>& operator=(SparseVector<T>&& rhs);

  /* dtor */
  virtual ~SparseVector() = default;

  /* iterator */
  iterator begin() { return a.begin(); }
  iterator end() { return a.end(); }

  iterator cbegin() const { return a.cbegin(); }
  iterator cend() const { return a.cend(); }

  /* I/O */
  template <typename U>
  friend std::ostream& operator<< (std::ostream& os, const SparseVector<U>& vin);

  std::ostream& print_compressed(std::ostream& os) const;

  /* indexing */
  reference operator[](size_type i);
  const_reference operator[] (size_type i) const;

  size_type get_ja(size_type j) const { return ja[j]; }

  /* dimension */
  size_type length() const { return n; }

  size_type nnonzero() const { return nnz; }

 protected:
  SparseVector(size_type n, size_type nnz): n(n), nnz(nnz), ja(nnz), a(nnz) {}

  size_type n;
  size_type nnz;

  FullVector<size_type> ja;
  FullVector<T> a;
};

IAJA_NAMESPACE_CLOSE

#include <iaja/linalg_vector.impl.h>

#endif //_LINALG_VECTOR_H_
