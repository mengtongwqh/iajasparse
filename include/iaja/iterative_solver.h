#include <iaja/linalg_matrix.h>
#include <iaja/linalg_vector.h>
#include <iaja/ilu.h>

#include <string>

class IterativeSolverILU {

 public:
  /* constructor and destructor */
  IterativeSolverILU(SparseMatrix<double>& A,
          unsigned int max_iter = 10000,
          double tol = 1.0e-6)
      : max_iter(max_iter), tol(tol), iter_count(0), A(A), ilu(A) {}
  virtual ~IterativeSolverILU() = default;
  unsigned int get_iter_count() { return iter_count; }
  double get_residual_norm() { return residual_norm; }

 protected:
  const unsigned int max_iter; // default = 10000
  const double tol; // default = 1.0e-6
  unsigned int iter_count;
  double residual_norm;
  SparseMatrix<double>& A;
  SparseILU ilu;

  void residual(const FullVector<double>& b,
          const FullVector<double>& x,
          FullVector<double>& res);

 public:
  void symbolic_factor(const std::string& reorder_method = "natural",
          unsigned int max_level_of_fill = 1);

  void test_ilu_factorization(
          const std::string& reorder_method,
          unsigned int max_level_of_fill,
          const FullVector<double>& rhs);
};



// ------------------------------------------
// PCG solver with ILU preconditioner
// ------------------------------------------
class PCG : public IterativeSolverILU {

 public:
  /* ctor and dtor */
  PCG(SparseMatrix<double>& A,
          unsigned int max_iter = 10000, double tol = 1e-6)
    : IterativeSolverILU(A, max_iter, tol) {}
  virtual ~PCG() = default;

  /* methods */
  int iterative_solve(const FullVector<double>& b, FullVector<double>& soln);

 private:
  unsigned int interval_residual_recompute = 50;
};


// ------------------------------------------
// Orthomin(k) solver with ILU preconditioner
// ------------------------------------------
class Orthomin : public IterativeSolverILU {

 public:
  /* ctor and dtor */
  Orthomin(SparseMatrix<double>& A, unsigned int k_orth = 5,
          unsigned int max_iter = 10000, double tol = 1e-6)
    : IterativeSolverILU(A, max_iter, tol), k_orth(k_orth) {}
  virtual ~Orthomin() = default;

  /* methods */
  int iterative_solve(const FullVector<double>& b, FullVector<double>& soln);

 private:
  const unsigned int k_orth;
};
