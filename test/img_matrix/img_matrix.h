#ifndef _IMG_MATRIX_H_
#define _IMG_MATRIX_H_

#include <iaja/iaja_config.h>
#include <iaja/linalg_vector.h>
#include <iaja/linalg_matrix.h>
#include <iaja/ilu.h>
#include <iaja/iterative_solver.h>

#include <string>

IAJA_NAMESPACE_OPEN

// -------------------------------------
// Image matrix test base class
// -------------------------------------
class ImgMatrixTest  {
 public:
  using size_type = std::size_t;
  ImgMatrixTest(size_type N);
  virtual ~ImgMatrixTest() = default;
  void img_lhs();
  void img_rhs();
  void set_diag_dominant_Mmatrix();

  // parameters
  size_type n;
  size_type N;

  // matrix structure
  SparseMatrixIaja<double> lhs;
  FullVector<double> rhs;
  FullVector<double> x;
  // PCG solver;

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

IAJA_NAMESPACE_CLOSE

#endif
