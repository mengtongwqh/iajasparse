#ifndef _ILU_H_
#define _ILU_H_

#include <iaja/iaja_config.h>
#include <iaja/linalg_matrix.h>
#include <iaja/linalg_vector.h>

#include <string>
#include <vector>

IAJA_NAMESPACE_OPEN

class SparseILURow : public SparseVector<double> {

  friend class SparseILU;

 public:
  // size_type is inherited from SparseVector
  using size_type = SparseVector<double>::size_type;

  // ctors and dtors
  SparseILURow():SparseVector(), level_of_fill() {}
  SparseILURow(size_t n, size_t nnz)
      : SparseVector<double>(n, nnz), level_of_fill(nnz), diag(0) {}
  SparseILURow(const SparseILURow& rhs);
  SparseILURow(SparseILURow&& rhs);
  virtual ~SparseILURow() = default;

  // interfaces
  SparseVector<double>::size_type diagonal() {return diag;}

 protected:
  FullVector<unsigned int> level_of_fill;
  size_type diag;
};


class SparseILU {

 public:
  using size_type = SparseILURow::size_type;
  /* constructors and destructors */
  explicit SparseILU(SparseMatrix<double>& mtrx);
  SparseILU(SparseILU&& rhs);
  virtual ~SparseILU() = default;

  /* operator overloading */
  friend std::ostream& operator<< (std::ostream& os, const SparseILU& ilu);

  /* public methods */
  void reorder(const std::string& reorder_method = "natural");
  void analyse(const unsigned int max_level_of_fill);
  void factor();
  void solve(const FullVector<double>& b, FullVector<double>& x) const;
  void print_level_of_fill(std::ostream& os) const;

 protected:
  SparseMatrix<double>& A;
  size_type  n;
  FullVector<size_type> order_new2old; // order[new_order] = old_order
  FullVector<size_type> order_old2new; // invord[old_order] = new_order
  std::vector<SparseILURow> rows;

 private:
  void merge_linked_list(const size_type rowid,
          const unsigned int max_level_of_fill,
          const size_type list_begin,
          std::vector<size_type>& row_linked_list,
          std::vector<unsigned int>& row_level_of_fill,
          size_type& nnzrow);
};

IAJA_NAMESPACE_CLOSE

#endif //_ILU_H_
