#ifndef _LINALG_VECTOR_H_
#define _LINALG_VECTOR_H_

#include <iaja/global_defs.h>
#include <iostream>
#include <vector>

IAJA_NAMESPACE_OPEN

template <typename T> class SparseMatrix;

/* A little more than the a naked pointer to a block of heap */
template<typename T>
class FullVector {

 public:
  // constructor and destructor
  FullVector();
  explicit FullVector(unsigned int n);
  FullVector(unsigned int n, const T* a);
  FullVector(const FullVector<T>& v);
  FullVector(FullVector<T>&& v);
  virtual ~FullVector();

  // operator overloading
  template <typename U>
  friend std::ostream& operator<< (std::ostream& os, const FullVector<U>& x);

  T& operator[] (unsigned int i);
  const T& operator[] (unsigned int i) const;
  T operator* (const FullVector<T>& x) const;
  FullVector<T>& operator= (const FullVector<T>& rhs);
  FullVector<T>& operator= (const std::vector<T>& rhs);
  operator T*() { return a; }

  // member functions
  unsigned int length() const {return n;}
  T norm_l2() const;
  void saxpy(T alpha, const FullVector<T>& x, const FullVector<T>& y);
  void add(const FullVector<T>& b, const FullVector<T>& c);
  void minus(const FullVector<T>& b, const FullVector<T>& c);
  void multiply(const SparseMatrix<T>& A, const FullVector<T>& x);
  void multiply(const T& s);

 protected:
  unsigned int n;
  T* a;
};


template <typename T>
class SparseVector {

 public:
  // constructors and destructors
  SparseVector();
  SparseVector(unsigned int n , unsigned int nnz);
  SparseVector(unsigned int n, unsigned int nnz, const unsigned int* ja, const T* a);
  SparseVector(const SparseVector<T>& rhs);
  SparseVector(SparseVector<T>&& rhs);
  virtual ~SparseVector() = default;

  /* operator overloading */
  template <typename U>
  friend std::ostream& operator<< (std::ostream& os, const SparseVector<U>& vin);
  T& operator[](unsigned int i);
  const T& operator[] (unsigned int i) const;

  /* public interfaces */
  unsigned int length() const { return n; }
  unsigned int nonzeros() const { return nnz; }

 protected:
  unsigned int n;
  unsigned int nnz;
  FullVector<unsigned int> ja;
  FullVector<T> a;
};

IAJA_NAMESPACE_CLOSE
#include <iaja/linalg_vector.impl.h>


#endif //_LINALG_VECTOR_H_
