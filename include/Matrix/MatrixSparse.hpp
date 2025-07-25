/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Basic/WarningMacro.hpp"
#include "LinearOp/ALinearOp.hpp"
#include "Matrix/AMatrix.hpp"
#include <Eigen/src/Core/Matrix.h>

#ifndef SWIG
DISABLE_WARNING_PUSH
DISABLE_WARNING_COND_EXPR_CONSTANT
DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE
DISABLE_WARNING_DECLARATION_HIDE_GLOBAL
#  include <Eigen/Sparse>
DISABLE_WARNING_POP
#endif

namespace gstlrn
{
class EOperator;
} // namespace gstlrn

namespace gstlrn
{
/**
 * Sparse Matrix
 *
 * Handle a sparse matrix that can be symmetrical, square or not.
 * Storage relies on Eigen3 Library.
 * Storage relies on Eigen3 Library.
 */
class GSTLEARN_EXPORT MatrixSparse: public AMatrix, public virtual ALinearOp
{
public:
  MatrixSparse(int nrow = 0, int ncol = 0, int ncolmax = -1);
  MatrixSparse(const MatrixSparse& m);
  MatrixSparse& operator=(const MatrixSparse& m);
  virtual ~MatrixSparse();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// Cloneable interface
  IMPLEMENT_CLONING(MatrixSparse)

  //// Interface to AStringable
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for ALinearOp

  int getSize() const override { return getNRows(); }

  /// Interface for AMatrix

  /*! Returns if the current matrix is Sparse */
  bool isSparse() const override { return true; }
  /*! Returns if the matrix belongs to the MatrixSparse class (avoids dynamic_cast) */
  bool isDense() const override { return false; }

  /*! Get the value from a matrix cell */
  double getValue(int row, int col) const override;
  /*! Set the value for a matrix cell */
  void setValue(int irow, int icol, double value) override;
  /*! Modifies the contents of a matrix cell */
  void updValue(int irow,
                int icol,
                const EOperator& oper,
                double value) override;
  /*! Set the contents of a Column */
  void setColumn(int icol, const VectorDouble& tab) override;
  /*! Set the contents of a Column to a constant */
  void setColumnToConstant(int icol, double value) override;
  /*! Set the contents of a Row */
  void setRow(int irow, const VectorDouble& tab) override;
  /*! Set the contents of a Row to a constant*/
  void setRowToConstant(int irow, double value) override;
  /*! Set the contents of the (main) Diagonal */
  void setDiagonal(const VectorDouble& tab) override;
  /*! Set the contents of the (main) Diagonal to a constant value */
  void setDiagonalToConstant(double value = 1.) override;
  /*! Add a value to each matrix component */
  void addScalar(double v) override;
  /*! Add value to matrix diagonal */
  void addScalarDiag(double v) override;
  /*! Multiply each matrix component by a value */
  void prodScalar(double v) override;
  /*! Set all the values of the matrix at once */
  void fill(double value) override;
  /*! Multiply the matrix row-wise */
  void multiplyRow(const VectorDouble& vec) override;
  /*! Multiply the matrix column-wise */
  void multiplyColumn(const VectorDouble& vec) override;
  /*! Divide the matrix row-wise */
  void divideRow(const VectorDouble& vec) override;
  /*! Divide the matrix column-wise */
  void divideColumn(const VectorDouble& vec) override;

  void resetFromValue(int nrows, int ncols, double value) override;
  void resetFromArray(int nrows, int ncols, const double* tab, bool byCol = true) override;
  void resetFromVD(int nrows, int ncols, const VectorDouble& tab, bool byCol = true) override;
  void resetFromVVD(const VectorVectorDouble& tab, bool byCol = true) override;

  /*! Transpose the matrix and return it as a copy*/
  MatrixSparse* transpose() const override;

  /*! Multiply matrix 'x' by matrix 'y' and store the result in 'this' */
  void prodMatMatInPlace(const AMatrix* x,
                         const AMatrix* y,
                         bool transposeX = false,
                         bool transposeY = false) override;
  /*! Perform 'this' = 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' */
  void prodNormMatMatInPlace(const AMatrix* a,
                             const AMatrix* m,
                             bool transpose = false) override;
  /*! Perform 'this' = 't(A)' %*% ['vec'] %*% 'A' or 'A' %*% ['vec'] %*% 't(A)' */
  void prodNormMatVecInPlace(const AMatrix* a,
                             const VectorDouble& vec,
                             bool transpose = false) override;
  /*! Perform 'this' = 't(A)' %*% 'A' or 'A' %*% 't(A)' */
  void prodNormMatInPlace(const AMatrix* a,
                          bool transpose = false) override;
  /*! Perform 'this' = 'val1' * 'mat1' + 'val2' * 'mat2' + 'val3' * 'mat3' */
  void linearCombination(double val1,
                         const AMatrix* mat1,
                         double val2         = 1.,
                         const AMatrix* mat2 = nullptr,
                         double val3         = 1.,
                         const AMatrix* mat3 = nullptr) override;
  /*! Add a matrix (multiplied by a constant) */
  void addMatInPlace(const AMatrix& y, double cx = 1., double cy = 1.) override;

  /*! Extract the contents of the matrix */
  NF_Triplet getMatrixToTriplet(int shiftRow = 0, int shiftCol = 0) const override;

  MatrixSparse* getRowAsMatrixSparse(int irow, double coeff = 1.) const;
  MatrixSparse* getColumnAsMatrixSparse(int icol, double coeff = 1.) const;

#ifndef SWIG
  int addVecInPlaceEigen(const Eigen::Map<const Eigen::VectorXd>& xm,
                         Eigen::Map<Eigen::VectorXd>& ym) const;
#endif

  // Static functions
  static MatrixSparse* create(const MatrixSparse* mat);
  static MatrixSparse* create(int nrow, int ncol);
  static MatrixSparse* createFromTriplet(const NF_Triplet& NF_T,
                                         int nrow    = 0,
                                         int ncol    = 0,
                                         int nrowmax = -1);
  static MatrixSparse* Identity(int nrow, double value = 1.);
  static MatrixSparse* addMatMat(const MatrixSparse* x,
                                 const MatrixSparse* y,
                                 double cx = 1.,
                                 double cy = 1.);
  static MatrixSparse* diagVec(const VectorDouble& vec);
  static MatrixSparse* diagConstant(int number, double value = 1.);
  static MatrixSparse* diagMat(MatrixSparse* A, int oper_choice);
  static MatrixSparse* glue(const MatrixSparse* A1,
                            const MatrixSparse* A2,
                            bool flagShiftRow,
                            bool flagShiftCol);
  static void glueInPlace(MatrixSparse* A1,
                          const MatrixSparse* A2,
                          bool flagShiftRow,
                          bool flagShiftCol);
  /*! Dump a specific range of samples from the internal storage */
  static void dumpElements(const String& title, int ifrom, int ito);

  void resetFromTriplet(const NF_Triplet& NF_T);

  /*! Set all the values of the Matrix with random values */
  void fillRandom(int seed = 432432, double zeroPercent = 0);

#ifndef SWIG
  int addVecInPlace(const constvect x, vect y) const;
#endif
  void addValue(int row, int col, double value);

  double L1Norm() const;
  void getStats(int* nrows, int* ncols, int* count, double* percent) const;
  int scaleByDiag();
  int addVecInPlaceVD(const VectorDouble& x, VectorDouble& y) const;
  void setConstant(double value);
  VectorDouble extractDiag(int oper_choice = 1) const;
  void prodNormDiagVecInPlace(const VectorDouble& vec, int oper = 1);
  MatrixSparse* extractSubmatrixByRanks(const VectorInt& rank_rows,
                                        const VectorInt& rank_cols) const;
  MatrixSparse* extractSubmatrixByColor(const VectorInt& colors,
                                        int ref_color,
                                        bool row_ok,
                                        bool col_ok);
  VectorInt colorCoding() const;
  int getNonZeros() const { return _getMatrixPhysicalSize(); }
  void gibbs(int iech, const VectorDouble& zcur, double* yk, double* sk);

  int forwardLU(const VectorDouble& b, VectorDouble& x, bool flagLower = true) const;
  void forceDimension(int maxRows, int maxCols);

#ifndef SWIG

protected:
  int _addToDest(const constvect inv, vect outv) const override;
#endif

#ifndef SWIG

public:
  void setDiagonal(const Eigen::Map<const Eigen::VectorXd>& tab);
  void setDiagonal(const constvect tab);
#endif

protected:
  void _allocate() override;
  void _deallocate() override;
  int _getIndexToRank(int irow, int icol) const override;
  void _setValueByRank(int rank, double value) override;
  void _transposeInPlace() override;
  int _invert() override;
  int _solve(const VectorDouble& b, VectorDouble& x) const override;
  int _getMatrixPhysicalSize() const override;
  double _getValueByRank(int rank) const override;
  double& _getValueRef(int irow, int icol) override;

#ifndef SWIG
  void _addProdVecMatInPlacePtr(constvect x, vect y, bool transpose = false) const override;
  void _addProdMatVecInPlacePtr(constvect x, vect y, bool transpose = false) const override;
#endif

  bool _isPhysicallyPresent(int /*irow*/, int /*icol*/) const override { return true; }
  void _setValues(const double* values, bool byCol) override;
  void _clear() override;
  bool _isElementPresent(int irow, int icol) const;
  void _allocate(int nrow, int ncol, int ncolmax);

private:
  static void _forbiddenForSparse(const String& func);
  int _eigen_findColor(int imesh,
                       int ncolor,
                       VectorInt& colors,
                       VectorInt& temp) const;

#ifndef SWIG

public:
  const Eigen::SparseMatrix<double>& eigenMat() const
  {
    return _eigenMatrix;
  }

  Eigen::SparseMatrix<double>& eigenMat()
  {
    return _eigenMatrix;
  }
#endif

private:
#ifndef SWIG
  Eigen::SparseMatrix<double> _eigenMatrix; // Storage (always stored Eigen::ColMajor)
#endif
  int _nColMax;
};

/*! Transform any matrix into a Sparse format */
GSTLEARN_EXPORT MatrixSparse* createFromAnyMatrix(const AMatrix* mat);

/*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' */
GSTLEARN_EXPORT MatrixSparse* prodNormMatMat(const MatrixSparse* a,
                                             const MatrixSparse* m,
                                             bool transpose = false);
/*! Product 't(A)' %*% 'vec' %*% 'A' or 'A' %*% 'vec' %*% 't(A)' stored in 'this'*/
GSTLEARN_EXPORT MatrixSparse* prodNormMatVec(const MatrixSparse* a,
                                             const VectorDouble& vec,
                                             bool transpose = false);
/*! Product 't(A)' %*% 'A' or 'A' %*% 't(A)' stored in 'this'*/
GSTLEARN_EXPORT MatrixSparse* prodNormMat(const MatrixSparse* a,
                                          bool transpose = false);
/*! Product 'Diag(vec)' %*% 'A' %*% 'Diag(vec)' */
GSTLEARN_EXPORT MatrixSparse* prodNormDiagVec(const MatrixSparse* a,
                                              const VectorDouble& vec,
                                              int oper_choice = 1);

// Not exported method

#ifndef SWIG
GSTLEARN_EXPORT Eigen::SparseMatrix<double> AtMA(const Eigen::SparseMatrix<double>& A,
                                                 const Eigen::SparseMatrix<double>& M);
#endif
} // namespace gstlrn