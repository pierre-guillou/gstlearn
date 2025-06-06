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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Matrix/AMatrix.hpp"
#include "Basic/WarningMacro.hpp"

#ifndef SWIG
DISABLE_WARNING_PUSH
DISABLE_WARNING_COND_EXPR_CONSTANT
DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE
DISABLE_WARNING_MAYBE_UNINITIALIZED
DISABLE_WARNING_DECLARATION_HIDE_GLOBAL
#  include <Eigen/Dense>
#  include <Eigen/Eigenvalues>
DISABLE_WARNING_POP
#endif

class MatrixSquare;
class MatrixSymmetric;
class EOperator;

/**
 * Dense Matrix
 * This class provides all the functions that can be performed using a Matrix stored
 * in "Dense" format (in opposition to the "Sparse" format).
 * This class can be derived in the case the matrix is Square, and even more if it is
 * Square and Symmetric.
 */

class GSTLEARN_EXPORT MatrixDense: public AMatrix
{

public:
  MatrixDense(int nrow = 0, int ncol = 0);
  MatrixDense(const MatrixDense& r);
  MatrixDense(const AMatrix& r);
  MatrixDense& operator=(const MatrixDense& r);
  virtual ~MatrixDense();

  /// Has a specific implementation in the Target language
  DECLARE_TOTL;

  /// Cloneable interface
  IMPLEMENT_CLONING(MatrixDense)

  /// Interface for AMatrix
  /*! Returns if the matrix belongs to the MatrixSparse class (avoids dynamic_cast) */
  bool isDense() const override { return true; }
  /*! Returns if the current matrix is Sparse */
  bool isSparse() const override { return false; }
  /*! Say if the matrix must be symmetric */
  bool mustBeSymmetric() const override { return false; }

  /*! Set the value for in a matrix cell */
  void setValue(int irow, int icol, double value, bool flagCheck = false) override;
  /*! Get the value from a matrix cell */
  virtual double getValue(int irow, int icol, bool flagCheck = false) const override;
  /*! Update the contents of a matrix cell */
  void updValue(int irow,
                int icol,
                const EOperator& oper,
                double value,
                bool flagCheck = false) override;

  /*! Set the contents of a Column */
  virtual void setColumn(int icol,
                         const VectorDouble& tab,
                         bool flagCheck = false) override;
  /*! Set the contents of a Column to a constant value */
  virtual void setColumnToConstant(int icol,
                                   double value,
                                   bool flagCheck = false) override;
  /*! Set the contents of a Row */
  virtual void setRow(int irow,
                      const VectorDouble& tab,
                      bool flagCheck = false) override;
  /*! Set the contents of a Row to a constant value*/
  virtual void setRowToConstant(int irow,
                                double value,
                                bool flagCheck = false) override;
  /*! Set the contents of the (main) Diagonal */
  virtual void setDiagonal(const VectorDouble& tab, bool flagCheck = false) override;
  /*! Set the contents of the (main) Diagonal to a constant value */
  virtual void setDiagonalToConstant(double value = 1.) override;
  /*! Add a value to each matrix component */
  virtual void addScalar(double v) override;
  /*! Add value to matrix diagonal */
  virtual void addScalarDiag(double v) override;
  /*! Multiply each matrix component by a value */
  virtual void prodScalar(double v) override;
  /*! Set all the values of the Matrix at once */
  virtual void fill(double value) override;
  /*! Multiply a Matrix row-wise */
  virtual void multiplyRow(const VectorDouble& vec) override;
  /*! Multiply a Matrix column-wise */
  virtual void multiplyColumn(const VectorDouble& vec) override;
  /*! Divide a Matrix row-wise */
  virtual void divideRow(const VectorDouble& vec) override;
  /*! Divide a Matrix column-wise */
  virtual void divideColumn(const VectorDouble& vec) override;
  /*! Perform 'vec' * 'this' */
  virtual VectorDouble prodVecMat(const VectorDouble& x, bool transpose = false) const override;
  /*! Perform 'this' * 'vec'*/
  virtual VectorDouble prodMatVec(const VectorDouble& x, bool transpose = false) const override;
  /*! Extract a Row */
  virtual VectorDouble getRow(int irow) const override;
  /*! Extract a Column */
  virtual VectorDouble getColumn(int icol) const override;
  constvect getColumnPtr(int icol) const;
  /*! Multiply matrix 'x' by matrix 'y' and store the result in 'this' */
  virtual void prodMatMatInPlace(const AMatrix* x,
                                 const AMatrix* y,
                                 bool transposeX = false,
                                 bool transposeY = false) override;

  static MatrixDense* create(const MatrixDense* mat);
  static MatrixDense* create(int nrow, int ncol);

  static MatrixDense* createFromVVD(const VectorVectorDouble& X);
  static MatrixDense* createFromVD(const VectorDouble& X,
                                   int nrow,
                                   int ncol,
                                   bool byCol             = false,
                                   bool invertColumnOrder = false);

  /// The next functions use specific definition of matrix (to avoid dynamic_cast)
  /// rather than manipulating AMatrix. They are not generic of AMatrix anymore.
  /// WARNING: output matrix should not match any of input matrices (speed up).

  /*! Add a matrix (multiplied by a constant) */
  void addMatInPlace(const MatrixDense& y, double cx = 1., double cy = 1.);
  /*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' stored in 'this'*/
  virtual void prodNormMatMatInPlace(const MatrixDense* a,
                                     const MatrixDense* m,
                                     bool transpose = false);
  /*! Product 't(A)' %*% ['vec'] %*% 'A' or 'A' %*% ['vec'] %*% 't(A)' stored in 'this'*/
  virtual void prodNormMatVecInPlace(const MatrixDense& a,
                                     const VectorDouble& vec = VectorDouble(),
                                     bool transpose          = false);

  VectorDouble getEigenValues() const { return _eigenValues; }
  const MatrixSquare* getEigenVectors() const { return _eigenVectors; }
  int invert2(MatrixDense& res) const;

  static MatrixDense* glue(const AMatrix* A1,
                           const AMatrix* A2,
                           bool flagShiftRow,
                           bool flagShiftCol);
  static bool sample(MatrixDense& res,
                     const AMatrix& A,
                     const VectorInt& rowKeep = VectorInt(),
                     const VectorInt& colKeep = VectorInt(),
                     bool flagInvertRow       = false,
                     bool flagInvertCol       = false);

  void unsample(const AMatrix* A,
                const VectorInt& rowFetch,
                const VectorInt& colFetch,
                bool flagInvertRow = false,
                bool flagInvertCol = false);

  MatrixDense compressMatLC(const MatrixDense& matLC, bool transpose = false);

  // Adding a Row or a Column (at the bottom or right of Rectangular Matrix)
  void addRow(int nrow_added = 1);
  void addColumn(int ncolumn_added = 1);

#ifndef SWIG
  static void sum(const MatrixDense* mat1,
                  const MatrixDense* mat2,
                  MatrixDense* mat3);
#endif
  void prodMatMatInPlaceOptim(const MatrixDense* x,
                              const MatrixDense* y,
                              bool transposeX = false,
                              bool transposeY = false);

protected:
  virtual void _allocate() override;
  virtual void _deallocate() override;

  virtual int _getMatrixPhysicalSize() const override;
  virtual double& _getValueRef(int irow, int icol) override;
  virtual double _getValueByRank(int rank) const override;
  virtual void _setValueByRank(int rank, double value) override;
  virtual int _getIndexToRank(int irow, int icol) const override;
  virtual void _transposeInPlace() override;
  virtual void _prodMatVecInPlacePtr(const double* x, double* y, bool transpose = false) const override;
  virtual void _prodVecMatInPlacePtr(const double* x, double* y, bool transpose = false) const override;
  virtual void _addProdMatVecInPlaceToDestPtr(const double* x, double* y, bool transpose = false) const override;
  virtual int _invert() override;
  virtual int _solve(const VectorDouble& b, VectorDouble& x) const override;

  int _computeEigen(bool optionPositive = true);
  int _computeGeneralizedEigen(const MatrixSymmetric& b, bool optionPositive = true);

private:
  void _recopy(const MatrixDense& r);
  int _terminateEigen(const Eigen::VectorXd& eigenValues,
                      const Eigen::MatrixXd& eigenVectors,
                      bool optionPositive = true,
                      bool changeOrder    = false);
  bool _needToReset(int nrows, int ncols) override;

#ifndef SWIG

public:
  constvect getViewOnColumn(int icol) const;
  vect getViewOnColumnModify(int icol);
#endif
#ifndef SWIG

public:
  Eigen::Map<const Eigen::MatrixXd> getEigenMat() const
  {
    return Eigen::Map<const Eigen::MatrixXd>(_eigenMatrix.data(), getNRows(), getNCols());
  }
  Eigen::Map<Eigen::MatrixXd> getEigenMat()
  {
    return Eigen::Map<Eigen::MatrixXd>(_eigenMatrix.data(), getNRows(), getNCols());
  }
#endif

protected:
  bool _flagEigenDecompose;
  VectorDouble _eigenValues;   // Used only when ! flag_eigen()
  MatrixSquare* _eigenVectors; // Used only when ! flag_eigen()
  int _maxSize;

protected:
  VectorDouble _eigenMatrix; // Eigen storage for Dense matrix in Eigen Library
};
