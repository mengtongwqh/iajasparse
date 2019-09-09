#ifndef _ILU_H_
#define _ILU_H_

#include <iaja/global_defs.h>
#include <iaja/linalg_matrix.h>
#include <iaja/linalg_vector.h>

#include <string>
#include <vector>

IAJA_NAMESPACE_OPEN

class SparseILURow : public SparseVector<double>{

  friend class SparseILU;

 public:
  SparseILURow():SparseVector(), level_of_fill() {}
  SparseILURow(unsigned int n, unsigned int nnz)
      : SparseVector<double>(n, nnz), level_of_fill(nnz), diag(0) {}
  SparseILURow(const SparseILURow& rhs);
  SparseILURow(SparseILURow&& rhs);
  virtual ~SparseILURow() = default;
  unsigned int diagonal() {return diag;}

 protected:
  FullVector<unsigned int> level_of_fill;
  unsigned int diag;
};


class SparseILU {

 public:
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
  unsigned int  n;
  FullVector<unsigned int> order_new2old; // order[new_order] = old_order
  FullVector<unsigned int> order_old2new; // invord[old_order] = new_order
  std::vector<SparseILURow> rows;

 private:
  void merge_linked_list(const unsigned int rowid,
          const unsigned int max_level_of_fill,
          const unsigned int list_begin,
          std::vector<unsigned int>& row_linked_list,
          std::vector<unsigned int>& row_level_of_fill,
          unsigned int& nnzrow);
};

IAJA_NAMESPACE_CLOSE

#endif //_ILU_H_
