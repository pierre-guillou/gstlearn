/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "Basic/Vector.hpp"
#include "Basic/IClonable.hpp"
#include "Basic/AStringable.hpp"
#include "Matrix/AMatrix.hpp"

/**
 * Rectangular matrices are stored by columns
 */
class GSTLEARN_EXPORT MatrixInt : public AStringable, public IClonable {

public:
  MatrixInt(int nrow = 0, int ncol = 0);
  MatrixInt(const MatrixInt &m);
  MatrixInt& operator= (const MatrixInt &r);
	virtual ~MatrixInt();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /*! Clonable interface */
  virtual IClonable* clone() const override { return new MatrixInt(*this); };

  void   reset(int nrows, int ncols);
  double getValue(int irow, int icol) const;
  double getValue(int irank) const;
  void   setValue(int rank, double value);
  void   setValue(int irow, int icol, double value);
  int    getMatrixSize() const;
  int    size() const { return getMatrixSize(); }
  VectorInt getValues() const;
  VectorVectorInt getMatrix() const;
  void   setValues(const VectorInt& values, bool byCol = true);
  void   setValues(const int* values, bool byCol = true);
  void   transposeInPlace();
  bool   empty() const { return _nRows <= 0 || _nCols <= 0; }

  int getNCols() const { return _nCols; }
  void setNCols(int cols) { _nCols = cols; }
  int getNRows() const { return _nRows; }
  void setNRows(int rows) { _nRows = rows; }

private:
  int    _getIndexToRank(int irow,int icol) const;
  void   _allocate();
  void   _deallocate();
  bool   _isIndexValid(int irow, int icol) const;
  bool   _isRankValid(int rank) const;
  bool   _isNumbersValid(int nrows, int ncols) const;

private:
  int _nRows;
  int _nCols;
  VectorInt _rectMatrix;
};