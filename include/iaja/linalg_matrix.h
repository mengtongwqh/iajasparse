#ifndef _LINALG_MATRIX_H_
#define _LINALG_MATRIX_H_

#include <iaja/iaja_config.h>
#include <iaja/linalg_vector.h>
#include <cstddef>
#include <iostream>

IAJA_NAMESPACE_OPEN

// ia-ja row compressed data structure
template <typename SizeType>
class SparsityPatternIaja {
 public:
  using size_type = SizeType;
  SparsityPatternIaja();
  ~SparsityPatternIaja();

  // operator overloading
  // friend std::ostream& operator<<(std::ostream& os, const SparsityPatternIaja& sp_pat);

  // matrix dimension
  size_type nrow() const { return nr; }
  size_type ncol() const { return nc; }
  size_type nnonzero() const { return nnz; }

  // display
  // std::ostream print_compressed();

 protected:
  FullVector<size_type> ia;
  FullVector<size_type> ja;
  size_type nr;
  size_type nc;
  size_type nnz;
};

/* sparse matrix in ia-ja (row compressed) data structure */
template <typename T>
class SparseMatrixIaja {
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
  SparseMatrixIaja(size_type nr, size_type nc, size_type nnz,
          const size_type* ia, const size_type *ja, const double *a);
  SparseMatrixIaja(size_type nr, size_type nc, size_type nnz);
  SparseMatrixIaja(const SparseMatrixIaja<T>& mat);
  virtual ~SparseMatrixIaja() = default;

  // operator overloading
  template <typename U>
  friend std::ostream& operator << (std::ostream& os, const SparseMatrixIaja<U>& mtrx);

  SparseMatrixIaja& operator=(const SparseMatrixIaja<T>& rhs);
  SparseMatrixIaja& operator=(const SparseMatrixIaja<T>&& rhs);

  T& operator[](size_type i);
  const T& operator[](size_type i) const;

  FullVector<T> operator*(const FullVector<T> x) const;

  // dimensions
  size_type nrow() const { return nr; }
  size_type ncol() const { return nc; }
  size_type nnonzero() const { return nnz; }

  SparseMatrixIaja<T> transpose() const;
  // void SparseMatrixIaja<T>::compress_storage();
  std::ostream& print_compressed(std::ostream& os) const;

  // ATTRIBUTES
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
