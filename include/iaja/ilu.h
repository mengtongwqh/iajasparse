#ifndef _ILU_H_
#define _ILU_H_

#include <iaja/iaja_config.h>
#include <iaja/linalg_matrix.h>
#include <iaja/linalg_vector.h>

#include <string>
#include <vector>

IAJA_NAMESPACE_OPEN

class IncompleteFactor {

 /* ================================== *
  * Abstract base class for all
  * incomplete factorization procedures
  * ================================== */

 public:
  class SparseFactorRow : public SparseVector<FloatType> {
   public:
    // size_type is inherited from SparseVector
    using size_type = SparseVector<FloatType>::size_type;
    // -----------------------------
    /* ctors */
    SparseFactorRow() : SparseVector(), level_of_fill() {}
    // SparseFactorRow(size_type n, size_type nnz):
        // SparseVector<FloatType>(n, nnz), level_of_fill(nnz), diag(0) {}
    SparseFactorRow(size_type n, FullVector<size_type>&& ja,
            FullVector<unsigned int>&& lof);
    /* copy ctor/assign */
    SparseFactorRow(const SparseFactorRow& rhs) = default;
    SparseFactorRow& operator= (const SparseFactorRow& rhs) = default;
    /* move ctor/assign */
    SparseFactorRow(SparseFactorRow&& rhs);
    SparseFactorRow& operator= (SparseFactorRow&& rhs);
    /* dtor */
    virtual ~SparseFactorRow() = default;
    // -----------------------------
    /* attributes */
    FullVector<unsigned int> level_of_fill;
    size_type diag;
  };

 public:
  using size_type = SparseFactorRow::size_type;
  /* ctor */
  IncompleteFactor(SparseMatrixIaja<FloatType>& mtrx,
          const std::string& reorder_method);
  /* copy ctor/assign */
  IncompleteFactor(const IncompleteFactor& rhs) = delete;
  IncompleteFactor& operator = (const IncompleteFactor& rhs) = delete;
  /* move ctor/assign */
  IncompleteFactor(IncompleteFactor&& rhs);
  IncompleteFactor& operator = (IncompleteFactor&& rhs);
  /* dtor */
  ~IncompleteFactor() = default;

  /* Numerical Operations */
  void reorder(const std::string& reorder_method = "natural");
  void analyse(const unsigned int max_Lof);
  virtual void factor();
  virtual void solve(const FullVector<FloatType>& b, FullVector<FloatType>& x) const;

  /* I/O */
  void print_level_of_fill(std::ostream& os) const;
  friend std::ostream& operator<< (std::ostream& os, const IncompleteFactor& ifac);

  /* Type Conversion */
  explicit operator SparseMatrixIaja<FloatType>() const;

 protected:
  SparseMatrixIaja<FloatType>& A;
  size_type n;
  FullVector<size_type> order_new2old;
  FullVector<size_type> order_old2new;
  std::vector<SparseFactorRow> rows;

 private:
  void merge_linked_list(const size_type rowid,
          const unsigned int max_level_of_fill,
          const size_type list_begin,
          std::vector<size_type>& row_linked_list,
          std::vector<unsigned int>& row_level_of_fill,
          size_type& nnzrow);
};



class SparseILU : public IncompleteFactor {

 /* ================================== *
  *  incomplete LU-factorization
  * ================================== */

 public:
  using size_type = SparseFactorRow::size_type;

  /* ctors */
  SparseILU(SparseMatrixIaja<FloatType>& mtrx, const std::string& reorder_method);
  /* copy ctor/assign */
  SparseILU(const SparseILU& rhs) = delete;
  SparseILU& operator = (const SparseILU&) = delete;
  /* move ctor/assign */
  SparseILU(SparseILU&& rhs);
  SparseILU& operator = (SparseILU&& rhs);
  /* dtor */
  virtual ~SparseILU() = default;

  /* Numeric Methods */
  virtual void factor();
  virtual void solve(const FullVector<FloatType>& b, FullVector<FloatType>& x) const;
};



class SparseIchol : public IncompleteFactor {

 /* ================================== *
  * incomplete Cholesky factorization
  * ================================== */

 public:
  /* ctor */
  SparseIchol(SparseMatrixIaja<FloatType>& mtrx, const std::string& reorder_method);
  /* copy ctor/assign */
  SparseIchol(const SparseIchol& rhs) = delete;
  SparseIchol& operator = (const SparseIchol& rhs) = delete;
  /* move ctor/assign */
  SparseIchol(SparseIchol&& rhs);
  SparseIchol& operator = (SparseIchol&& rhs);
  /* dtor */
  virtual ~SparseIchol() = default;

  /* Numerical Operations */
  virtual void factor();
  virtual void solve(const FullVector<FloatType>& b, FullVector<FloatType>& x) const;
};

IAJA_NAMESPACE_CLOSE

#endif //_ILU_H_
