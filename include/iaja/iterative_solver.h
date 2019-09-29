#include <iaja/iaja_config.h>
#include <iaja/linalg_matrix.h>
#include <iaja/linalg_vector.h>
#include <iaja/incomplete_factor.h>

#include <string>

IAJA_NAMESPACE_OPEN

class IterativeSolver {
 public:
  /* ctor */
  IterativeSolver(SparseMatrixIaja<FloatType>& A,
          unsigned int max_iter = 10000,
          FloatType tol = 1.0e-6):
      A(A), max_iter(max_iter), tol(tol),
      iter_count(0), residual_norm(0.0) {};
  /* dtor */
  ~IterativeSolver() = default;
    
  /* solver state accessor */
  unsigned int get_iter_count() { return iter_count; }
  double get_residual_norm() { return residual_norm; }
  virtual int iterative_solve(const FullVector<FloatType>& b, FullVector<FloatType>& x);

 protected:
  void residual(const FullVector<double>& b,
          const FullVector<double>& x,
          FullVector<double>& res);

  SparseMatrixIaja<FloatType>& A;
  unsigned int max_iter;
  unsigned int iter_count;
  FloatType tol;
  FloatType residual_norm;
};


class IterativeSolverIFactor : public IterativeSolver {

 public:
  /* ctor */
  IterativeSolverIFactor(IncompleteFactor* ifac,
          unsigned int max_iter = 10000,
          double tol = 1.0e-6):
    IterativeSolver(ifac->A, max_iter, tol), ifac(ifac) {}
  /* dtor */
  virtual ~IterativeSolverIFactor() = default;

  /* accessor */
  const SparseILU& get_ilu() const { return *ilu; }

  /* numerical operations */
  virtual void symbolic_factor(const std::string& reorder_method = "natural",
          unsigned int max_level_of_fill = 1);

  void test_incomplete_factorization(
          const std::string& reorder_method,
          unsigned int max_level_of_fill,
          const FullVector<double>& rhs);

 protected:
  IncompleteFactor* ifac;
};



// ------------------------------------------
// PCG solver with ILU preconditioner
// ------------------------------------------
class PCG : public IterativeSolverIFactor {

 public:
  /* ctor and dtor */
  PCG(IncompleteFactor* ifac,
          unsigned int max_iter = 10000, double tol = 1e-6)
    : IterativeSolverIFactor(ifac, max_iter, tol) {}
  virtual ~PCG() = default;

  /* methods */
  int iterative_solve(const FullVector<double>& b, FullVector<double>& soln);

 private:
  const unsigned int interval_residual_recompute = 5;
};


// ------------------------------------------
// Orthomin(k) solver with ILU preconditioner
// ------------------------------------------
class Orthomin : public IterativeSolverIFactor {

 public:
  /* ctor and dtor */
  Orthomin(IncompleteFactor* ifac, unsigned int k_orth = 5,
          unsigned int max_iter = 10000, double tol = 1e-6)
    : IterativeSolverIFactor(ifac, max_iter, tol), k_orth(k_orth) {}
  virtual ~Orthomin() = default;

  /* methods */
  int iterative_solve(const FullVector<double>& b, FullVector<double>& soln);

 private:
  const unsigned int k_orth;
};

IAJA_NAMESPACE_CLOSE
