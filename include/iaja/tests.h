#ifndef _IMG_MATRIX_H_
#define _IMG_MATRIX_H_

#include <iaja/iaja_config.h>
#include <iaja/linalg_vector.h>
#include <iaja/linalg_matrix.h>
#include <iaja/ilu.h>
#include <iaja/iterative_solver.h>

#include <string>

IAJA_NAMESPACE_OPEN
/* ==================================================== *
 *            IMAGE DENOISING TEST PROBLEM              *
 * ==================================================== */

// -------------------------------------
// Image matrix test base class
// -------------------------------------
class ImgMatrixTest  {
 public:
  using size_type = SizeType;
  explicit ImgMatrixTest(size_type N);
  virtual ~ImgMatrixTest() = default;

  void img_lhs();
  void img_rhs();
  void set_diag_dominant_Mmatrix();

  // parameters
  size_type n;
  size_type N;

  // matrix structure
  SparseMatrixIaja<FloatType> lhs;
  FullVector<FloatType> rhs;
  FullVector<FloatType> x;

 protected:
  double alpha;
};

// -------------------------------------
// Image matrix test w. PCG
// -------------------------------------
class ImgMatrixPCG : public ImgMatrixTest {

 public:
  ImgMatrixPCG(size_type N, unsigned int max_iter, double tol)
      : ImgMatrixTest(N), solver(lhs, max_iter, tol) {}
  virtual ~ImgMatrixPCG() = default;
  void test_ilu_procedures(unsigned int max_level_of_fill, const std::string& reorder_method);

  PCG solver;
};

// -------------------------------------
// Image matrix test w. Orthomin
// -------------------------------------
class ImgMatrixOrthomin : public ImgMatrixTest {
 public:
  ImgMatrixOrthomin(size_type N, unsigned int k_orth,
          unsigned int max_iter, double tol):
      ImgMatrixTest(N), solver(lhs, k_orth, max_iter, tol) {}
  virtual ~ImgMatrixOrthomin() = default;

  Orthomin solver;
};
/* ==================================================== */


/* ==================================================== *
 *           3D FINITE DIFFERENCE GRID TEST             *
 * ==================================================== */
class FDGrid {

 public:
  using size_type = SizeType;
  /* Ctor */
  explicit FDGrid(size_type n_dim, const std::string& method = "ascend");
  /* Copy Ctor/Assign */
  FDGrid(const FDGrid& rhs) = delete;
  FDGrid& operator = (const FDGrid& rhs) = delete;
  /* Move Ctor/Assign */
  FDGrid(FDGrid&& rhs) = delete;
  FDGrid& operator = (FDGrid&& rhs) = delete;
  /* Dtor */
  ~FDGrid() = default;

  /* getters */
  const SparseMatrixIaja<FloatType>& get_lhs() const {return lhs;}
  const FullVector<FloatType>& get_rhs() const {return rhs;}

 protected:
  size_type nx, ny, nz, n;
  FullVector<FloatType> rhs;
  FullVector<FloatType> x;
  SparseMatrixIaja<FloatType> lhs;
  FullVector<SizeType> idiag;

 private:
  void build_sparsity(FullVector<SizeType>& ia, FullVector<SizeType>& ja);
  void set_ascend();
  void set_anisotropic_K();
};

class FDGridPCG : public FDGrid, public PCG {
 public:
  FDGridPCG(size_type n_dim, const std::string& method = "ascend",
          unsigned int max_iter = 1000, double tol = 1e-6):
      FDGrid(n_dim, method), PCG(lhs, max_iter, tol) {}
  ~FDGridPCG() = default;

  int iterative_solve() { return PCG::iterative_solve(rhs, x); }
};


class FDGridOrthomin : public FDGrid, public Orthomin {
 public:
  FDGridOrthomin(size_type n_dim, const std::string& method = "ascend",
          unsigned int k_orth = 5, unsigned int max_iter = 1000, double tol = 1e-6):
      FDGrid(n_dim, method), Orthomin(lhs, k_orth, max_iter, tol) {}
  ~FDGridOrthomin() = default;

  int iterative_solve() { return Orthomin::iterative_solve(rhs, x); }
};


/* ==================================================== */
IAJA_NAMESPACE_CLOSE

#endif
