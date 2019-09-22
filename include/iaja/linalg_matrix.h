#ifndef _LINALG_MATRIX_H_
#define _LINALG_MATRIX_H_

#include <iaja/iaja_config.h>
#include <iaja/linalg_vector.h>
#include <iostream>

IAJA_NAMESPACE_OPEN

// ia-ja row compressed data structure
class SparsityIaja {

 public:
  using size_type = SizeType;

  /* ctor */
  SparsityIaja();
  // SparsityIaja(size_type nr, size_type nc, size_type nnz);
  SparsityIaja(size_type nr, size_type nc, size_type nnz, const size_type* ia, const size_type* ja);
  SparsityIaja(size_type nc, const FullVector<size_type>& ia, const FullVector<size_type>& ja);
  /* move ctor */
  SparsityIaja(SparsityIaja&& rhs);
  /* move assign */
  SparsityIaja& operator= (SparsityIaja&& rhs);
  /* copy ctor */
  SparsityIaja(const SparsityIaja& rhs) = default;
  /* copy assign */
  SparsityIaja& operator= (const SparsityIaja& rhs) = default;
  /* dtor */
  virtual ~SparsityIaja() = default;

  /* I/O */
  friend std::ostream& operator<<(std::ostream& os, const SparsityIaja& sp_pat);
  std::ostream& print_compressed(std::ostream& os);

  /* matrix dimension */
  size_type nrow() const { return nr; }
  size_type ncol() const { return nc; }
  size_type nnonzero() const { return nnz; }

 protected:
  size_type nr;
  size_type nc;
  size_type nnz;
  FullVector<size_type> ia;
  FullVector<size_type> ja;
};

/* sparse matrix in ia-ja (row compressed) data structure */
template <typename T>
class SparseMatrixIaja : public SparsityIaja {

  /* classes having direct access to sparsity structure */
  friend class SparseILU;

 public:
  using value_type      = T;
  using size_type       = SparsityIaja::size_type;
  using pointer         = value_type *;
  using const_pointer   = const value_type *;
  using iterator        = value_type *;
  using const_iterator  = const value_type *;
  using reference       = value_type &;
  using const_reference = const value_type &;

  /* ctor */
  SparseMatrixIaja() : SparsityIaja(), a() {}
  SparseMatrixIaja(size_type nr, size_type nc, size_type nnz,
          const size_type* ia, const size_type *ja);
  SparseMatrixIaja(size_type nr, size_type nc, size_type nnz,
          const size_type* ia, const size_type *ja, const T& a);
  SparseMatrixIaja(size_type nr, size_type nc, size_type nnz,
          const size_type* ia, const size_type *ja, const T* a);
  SparseMatrixIaja(size_type nc,
          const FullVector<size_type>& ia,
          const FullVector<size_type>& ja,
          const FullVector<T>& a);
  /* copy ctor */
  SparseMatrixIaja(const SparseMatrixIaja<T>& mat) = default;
  /* copy assign */
  SparseMatrixIaja<T>& operator= (const SparseMatrixIaja<T>& rhs) = default;
  /* move ctor */
  SparseMatrixIaja(SparseMatrixIaja<T>&& mat);
  /* move assign */
  SparseMatrixIaja<T>& operator=(SparseMatrixIaja<T>&& mat);
  /* dtor */
  virtual ~SparseMatrixIaja() = default;

  /* I/O */
  template <typename U>
  friend std::ostream& operator << (std::ostream& os, const SparseMatrixIaja<U>& mtrx);

  /* Indexing */
  reference operator[](size_type i);
  const_reference operator[](size_type i) const;

  FullVector<T> operator*(const FullVector<T> x) const;

  // dimensions
  size_type nrow() const { return nr; }
  size_type ncol() const { return nc; }
  size_type nnonzero() const { return nnz; }

  SparseMatrixIaja<T> transpose() const;
  // void SparseMatrixIaja<T>::compress_storage();
  std::ostream& print_compressed(std::ostream& os) const;

 protected:
  FullVector<T> a;        // nonzero entries
};

IAJA_NAMESPACE_CLOSE

#include <iaja/linalg_matrix.impl.h>

#endif  //_LINALG_MATRIX_H_
