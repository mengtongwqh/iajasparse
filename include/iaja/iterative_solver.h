#ifndef _ITERATIVE_SOLVER_H_
#define _ITERATIVE_SOLVER_H_

#include <iaja/iaja_config.h>
#include <iaja/linalg_matrix.h>
#include <iaja/linalg_vector.h>
#include <iaja/incomplete_factor.h>

#include <string>

IAJA_NAMESPACE_OPEN

class IterativeSolver {
 public:
  /* ctor */
  IterativeSolver(const SparseMatrixIaja<FloatType>& A,
          unsigned int max_iter = 1000,
          FloatType tol = 1.0e-6):
      A(A), max_iter(max_iter), iter_count(0),
      tol(tol), residual_norm(0.0) {}

  /* dtor */
  virtual ~IterativeSolver() = default;

  /* solver state accessor */
  unsigned int get_iter_count() { return iter_count; }
  double get_residual_norm() { return residual_norm; }
  virtual int iterative_solve(const FullVector<FloatType>& b, FullVector<FloatType>& x) = 0;

 protected:
  void residual(const FullVector<double>& b,
          const FullVector<double>& x,
          FullVector<double>& res);

  const SparseMatrixIaja<FloatType>& A;
  unsigned int max_iter;
  unsigned int iter_count;
  FloatType tol;
  FloatType residual_norm;
};


class IterativeSolverIFactor : public IterativeSolver {

 public:
  /* ctor */
  IterativeSolverIFactor(IncompleteFactor* ifac,
          unsigned int max_iter = 1000,
          double tol = 1.0e-6):
      IterativeSolver(ifac->get_A(), max_iter, tol), ifac(ifac) {}

  IterativeSolverIFactor(const SparseMatrixIaja<FloatType>& A,
          unsigned int max_iter = 1000,
          double tol = 1.0e-6):
      IterativeSolver(A, max_iter, tol), ifac(nullptr) {}

  /* dtor */
  virtual ~IterativeSolverIFactor() = default;

  /* accessor */
  const IncompleteFactor& get_ilu() const { return *ifac; }

  /* numerical operations */
  virtual void symbolic_factor(unsigned int level_of_fill = 1);

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
          unsigned int max_iter = 1000,
          double tol = 1.0e-6):
      IterativeSolverIFactor(ifac, max_iter, tol) {}

  PCG(const SparseMatrixIaja<FloatType>& A,
          unsigned int max_iter = 1000,
          double tol = 1.0e-6):
      IterativeSolverIFactor(A, max_iter, tol) {}

  virtual ~PCG() = default;

  /* methods */
  virtual int iterative_solve(const FullVector<double>& b, FullVector<double>& soln);

 private:
  const unsigned int interval_residual_recompute = 50;
};


// ------------------------------------------
// Orthomin(k) solver with ILU preconditioner
// ------------------------------------------
class Orthomin : public IterativeSolverIFactor {

 public:
  /* ctor */
  Orthomin(IncompleteFactor* ifac,
          unsigned int max_iter = 1000,
          double tol = 1.0e-6,
          unsigned int k_orth = 5):
      IterativeSolverIFactor(ifac, max_iter, tol), k_orth(k_orth) {}

  Orthomin(const SparseMatrixIaja<FloatType>& A,
          unsigned int max_iter = 1000,
          double tol = 1.0e-6,
          unsigned int k_orth = 5):
      IterativeSolverIFactor(A, max_iter, tol), k_orth(k_orth) {}

  /* dtor */
  virtual ~Orthomin() = default;

  /* methods */
  virtual int iterative_solve(const FullVector<double>& b, FullVector<double>& soln);

 private:
  const unsigned int k_orth;
};

IAJA_NAMESPACE_CLOSE

#endif
