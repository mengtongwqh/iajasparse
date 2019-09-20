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
  // standard data types similar to std library containers
  using value_type      = T;
  using size_type       = SizeType;
  using pointer         = value_type *;
  using const_pointer   = const value_type *;
  using iterator        = value_type *;
  using const_iterator  = const value_type *;
  using reference       = value_type &;
  using const_reference = const value_type &;

  // constructor and destructor
  FullVector();
  explicit FullVector(size_type n);
  FullVector(size_type n, const T* a);
  FullVector(const FullVector<T>& v);
  FullVector(FullVector<T>&& v);
  virtual ~FullVector();

  // operator overloading
  template <typename U>
  friend std::ostream& operator<< (std::ostream& os, const FullVector<U>& x);

  T& operator[] (size_type i);
  const T& operator[] (size_type i) const;
  T& operator[] (long int i);
  const T& operator[] (long int i) const;
  T& operator[] (int i);
  const T& operator[] (int i) const;

  T operator* (const FullVector<T>& x) const;
  FullVector<T>& operator= (const FullVector<T>& rhs);
  FullVector<T>& operator= (const std::vector<T>& rhs);
  operator T*() { return a; }

  // member functions
  size_type length() const {return n;}
  T norm_l2() const;
  void saxpy(T alpha, const FullVector<T>& x, const FullVector<T>& y);
  void add(const FullVector<T>& b, const FullVector<T>& c);
  void minus(const FullVector<T>& b, const FullVector<T>& c);
  void multiply(const SparseMatrixIaja<T>& A, const FullVector<T>& x);
  void multiply(const T& s);
  void cumsum();

 protected:
  size_type n;
  T* a;
};


template <typename T>
class SparseVector {

 public:
  // standard data types similar to std library containers
  using value_type      = T;
  using size_type       = SizeType;
  using pointer         = value_type *;
  using const_pointer   = const value_type *;
  using iterator        = value_type *;
  using const_iterator  = const value_type *;
  using reference       = value_type &;
  using const_reference = const value_type &;

  // constructors and destructors
  SparseVector();
  SparseVector(size_type n , size_type nnz);
  SparseVector(size_type n, size_type nnz, const size_type* ja, const T* a);
  SparseVector(const SparseVector<T>& rhs);
  SparseVector(SparseVector<T>&& rhs);
  virtual ~SparseVector() = default;

  /* operator overloading */
  template <typename U>
  friend std::ostream& operator<< (std::ostream& os, const SparseVector<U>& vin);
  T& operator[](size_type i);
  const T& operator[] (size_type i) const;

  /* public interfaces */
  size_type length() const { return n; }
  size_type nonzeros() const { return nnz; }

 protected:
  size_type n;
  size_type nnz;
  FullVector<size_type> ja;
  FullVector<T> a;
};

IAJA_NAMESPACE_CLOSE

#include <iaja/linalg_vector.impl.h>

#endif //_LINALG_VECTOR_H_
