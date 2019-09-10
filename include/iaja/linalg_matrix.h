#ifndef _LINALG_MATRIX_H_
#define _LINALG_MATRIX_H_

#include <iaja/iaja_config.h>
#include <iaja/linalg_vector.h>
#include <cstddef>
#include <iostream>

IAJA_NAMESPACE_OPEN

/* sparse matrix in ia-ja (row compressed) data structure */
template <typename T>
class SparseMatrix {

 public:
  // standard data types similar to std library containers
  using value_type      = T;
  using size_type       = std::size_t;
  using pointer         = value_type *;
  using const_pointer   = const value_type *;
  using iterator        = value_type *;
  using const_iterator  = const value_type *;
  using reference       = value_type &;
  using const_reference = const value_type &;

  // constructors and destructors
  SparseMatrix(size_type nr, size_type nc, size_type nnz,
          const size_type* ia, const size_type *ja, const double *a);
  SparseMatrix(size_type nr, size_type nc, size_type nnz);
  SparseMatrix(const SparseMatrix<T>& mat);
  virtual ~SparseMatrix() = default;

  // operator overloading
  template <typename U>
  friend std::ostream& operator << (std::ostream& os, const SparseMatrix<U>& mtrx);
  SparseMatrix& operator=(const SparseMatrix<T>& rhs);
  SparseMatrix& operator=(const SparseMatrix<T>&& rhs);
  T& operator[](size_type i);
  const T& operator[](size_type i) const;
  FullVector<T> operator*(const FullVector<T> x) const;

  // methods
  size_type nrow() const { return nr; }
  size_type ncol() const { return nc; }
  size_type nonzeros() const { return nnz; }
  SparseMatrix<T> transpose() const;

  FullVector<size_type> ia; // row index
  FullVector<size_type> ja; // column index
  FullVector<T> a;        // nonzero entries

 protected:
  size_type nr;
  size_type nc;
  size_type nnz; // count of nonzeros
};

IAJA_NAMESPACE_CLOSE

#include <iaja/linalg_matrix.impl.h>

#endif  //_LINALG_MATRIX_H_
