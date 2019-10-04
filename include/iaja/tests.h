#ifndef _TESTS_H_
#define _TESTS_H_

#include <iaja/iaja_config.h>
#include <iaja/incomplete_factor.h>
#include <iaja/iterative_solver.h>
#include <iaja/linalg_vector.h>
#include <iaja/linalg_matrix.h>

#include <string>

IAJA_NAMESPACE_OPEN

/* ==================================================== *
 *            IMAGE DENOISING TEST PROBLEM              *
 * ==================================================== */

// -------------------------------------
// Image matrix test base class
// -------------------------------------
class ImgDenoiseTest  {
 public:
  using size_type = SizeType;

  explicit ImgDenoiseTest(size_type N);
  virtual ~ImgDenoiseTest() = default;

  void img_lhs();
  void img_rhs();
  void set_diag_dominant_Mmatrix();

  const SparseMatrixIaja<FloatType>& get_lhs() {return lhs;}
  const FullVector<FloatType>& get_rhs() {return rhs;}

 protected:
  // parameters
  size_type n;
  size_type N;

  // matrix structure
  SparseMatrixIaja<FloatType> lhs;
  FullVector<FloatType> rhs;

 public:
  FullVector<FloatType> x;

 protected:
  double alpha;
};



/* ==================================================== *
 *           3D FINITE DIFFERENCE GRID TEST             *
 * ==================================================== */
class EllipticalFDTest {

 public:
  using size_type = SizeType;

  /* Ctor */
  explicit EllipticalFDTest(size_type n_dim, const std::string& method = "ascend");

  /* Copy Ctor/Assign */
  EllipticalFDTest(const EllipticalFDTest& rhs) = delete;
  EllipticalFDTest& operator = (const EllipticalFDTest& rhs) = delete;

  /* Move Ctor/Assign */
  EllipticalFDTest(EllipticalFDTest&& rhs) = delete;
  EllipticalFDTest& operator = (EllipticalFDTest&& rhs) = delete;

  /* Dtor */
  ~EllipticalFDTest() = default;

  /* getters */
  const SparseMatrixIaja<FloatType>& get_lhs() const {return lhs;}
  const FullVector<FloatType>& get_rhs() const {return rhs;}

 protected:
  size_type nx, ny, nz, n;
  FullVector<FloatType> rhs;
  SparseMatrixIaja<FloatType> lhs;
  FullVector<SizeType> idiag;

 public:
  FullVector<FloatType> x;

 private:
  void build_sparsity(FullVector<SizeType>& ia, FullVector<SizeType>& ja);
  void set_ascend();
  void set_anisotropic_K();
  FloatType standard_solution_ascend() const;
};

/* ==================================================== */

IAJA_NAMESPACE_CLOSE

#endif
