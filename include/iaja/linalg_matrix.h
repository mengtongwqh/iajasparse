#ifndef _LINALG_MATRIX_H_
#define _LINALG_MATRIX_H_

#include <iostream>
#include <iaja/linalg_vector.h>


/* sparse matrix in ia-ja (row compressed) data structure */
template <typename T>
class SparseMatrix {

 public:
  // constructors and destructors
  SparseMatrix(unsigned int nr, unsigned int nc, unsigned int nnz,
          const unsigned int* ia, const unsigned int *ja, const double *a);
  SparseMatrix(unsigned int nr, unsigned int nc, unsigned int nnz);
  SparseMatrix(const SparseMatrix<T>& mat);
  virtual ~SparseMatrix() = default;

  // operator overloading
  template <typename U>
  friend std::ostream& operator << (std::ostream& os, const SparseMatrix<U>& mtrx);
  SparseMatrix& operator=(const SparseMatrix<T>& rhs);
  SparseMatrix& operator=(const SparseMatrix<T>&& rhs);
  T& operator[](unsigned int i);
  const T& operator[](unsigned int i) const;
  FullVector<T> operator*(const FullVector<T> x) const;

  // methods
  unsigned int nrow() const { return nr; }
  unsigned int ncol() const { return nc; }
  unsigned int nonzeros() const { return nnz; }

  FullVector<unsigned int> ia; // row index
  FullVector<unsigned int> ja; // column index
  FullVector<T> a;        // nonzero entries

 protected:
  unsigned int nr;
  unsigned int nc;
  unsigned int nnz; // count of nonzeros
};

#include <iaja/linalg_matrix.impl.h>

#endif  //_LINALG_MATRIX_H_
