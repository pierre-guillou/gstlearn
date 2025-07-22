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
#include "Matrix/MatrixSparse.hpp"
#include "Basic/AException.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/Law.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/WarningMacro.hpp"
#include "Matrix/AMatrix.hpp"
#include "Matrix/LinkMatrixSparse.hpp"
#include "Matrix/MatrixFactory.hpp"
#include "Matrix/NF_Triplet.hpp"

#include <Eigen/src/Core/Map.h>
#include <Eigen/src/Core/Matrix.h>
#include <iostream>

DISABLE_WARNING_PUSH
DISABLE_WARNING_COND_EXPR_CONSTANT
DISABLE_WARNING_UNUSED_BUT_SET_VARIABLE
#include <Eigen/SparseCholesky>
#include <omp.h>
DISABLE_WARNING_POP

/**
 * This variable switches ON/OFF the ability to use Eigen library for Algebra
 */

namespace gstlrn
{

MatrixSparse::MatrixSparse(int nrow, int ncol, int ncolmax)
  : AMatrix(nrow, ncol)
  , _eigenMatrix()
  , _nColMax(ncolmax)
{
  _allocate();
}

MatrixSparse::MatrixSparse(const MatrixSparse& m)
  : AMatrix(m)
  // , ALinearOp(m)
  , _eigenMatrix(m._eigenMatrix)
{
}

MatrixSparse& MatrixSparse::operator=(const MatrixSparse& m)
{
  if (this != &m)
  {
    AMatrix::operator=(m);
    // ALinearOp::operator=(m);
    if (!m.empty())
    {
      _eigenMatrix = m._eigenMatrix;
    }
  }
  return *this;
}

MatrixSparse::~MatrixSparse()
{
  _deallocate();
}

void MatrixSparse::resetFromValue(int nrows, int ncols, double value)
{
  DECLARE_UNUSED(nrows);
  DECLARE_UNUSED(ncols);
  DECLARE_UNUSED(value);
  _forbiddenForSparse("resetFromValue");
}

void MatrixSparse::resetFromArray(int nrows, int ncols, const double* tab, bool byCol)
{
  DECLARE_UNUSED(nrows);
  DECLARE_UNUSED(ncols);
  DECLARE_UNUSED(tab);
  DECLARE_UNUSED(byCol);
  _forbiddenForSparse("resetFromArray");
}

void MatrixSparse::resetFromVD(int nrows, int ncols, const VectorDouble& tab, bool byCol)
{
  DECLARE_UNUSED(nrows);
  DECLARE_UNUSED(ncols);
  DECLARE_UNUSED(tab);
  DECLARE_UNUSED(byCol);
  _forbiddenForSparse("resetFromVD");
}

void MatrixSparse::resetFromVVD(const VectorVectorDouble& tab, bool byCol)
{
  DECLARE_UNUSED(tab);
  DECLARE_UNUSED(byCol);
  _forbiddenForSparse("resetFromVVD");
}

void MatrixSparse::resetFromTriplet(const NF_Triplet& NF_T)
{
  eigenMat() = NF_T.buildEigenFromTriplet();
  _setNRows(eigenMat().rows());
  _setNCols(eigenMat().cols());
}

void MatrixSparse::fillRandom(int seed, double zeroPercent)
{
  law_set_random_seed(seed);

  int nrow = getNRows();
  int ncol = getNCols();
  NF_Triplet NF_T;
  for (int irow = 0; irow < nrow; irow++)
    for (int icol = 0; icol < ncol; icol++)
    {
      if (law_uniform(0., 1.) < zeroPercent) continue;
      NF_T.add(irow, icol, law_gaussian());
    }
  NF_T.force(nrow, ncol);
  resetFromTriplet(NF_T);
}

void MatrixSparse::_transposeInPlace()
{
  Eigen::SparseMatrix<double> temp;
  temp = eigenMat().transpose();
  eigenMat().swap(temp);
}

MatrixSparse* MatrixSparse::transpose() const
{
  MatrixSparse* mat = dynamic_cast<MatrixSparse*>(clone());
  mat->eigenMat()   = eigenMat().transpose();
  return mat;
}

/**
 * Fill a column of an already existing Sparse matrix, using 'tab' as entry
 * The input 'tab' corresponds to the whole column contents
 * @param icol Column rank
 * @param tab  Vector containing the information (Dimension: nrows)
 */
void MatrixSparse::setColumn(int icol, const VectorDouble& tab)
{
  int nrows = getNRows();
  if (getFlagMatrixCheck())
  {
    if (!_isColumnValid(icol)) return;
    if (!_isColumnSizeConsistent(tab)) return;
  }
  for (int irow = 0; irow < nrows; irow++)
    eigenMat().coeffRef(irow, icol) = tab[irow];
}

void MatrixSparse::setColumnToConstant(int icol, double value)
{
  int nrows = getNRows();
  if (getFlagMatrixCheck())
  {
    if (!_isColumnValid(icol)) return;
  }
  for (int irow = 0; irow < nrows; irow++)
    eigenMat().coeffRef(irow, icol) = value;
}

/**
 * Fill a row of an already existing Sparse matrix, using 'tab' as entry
 * The input 'tab' corresponds to the whole row contents
 * @param irow Row rank
 * @param tab  Vector containing the information (Dimension: ncols)
 *
 * @warning: This method only copies the values at the non-zero existing entries
 */
void MatrixSparse::setRow(int irow, const VectorDouble& tab)
{
  int ncols = getNCols();
  if (getFlagMatrixCheck())
  {
    if (!_isRowValid(irow)) return;
    if (!_isRowSizeConsistent(tab)) return;
  }
  for (int icol = 0; icol < ncols; icol++)
    eigenMat().coeffRef(irow, icol) = tab[icol];
}

void MatrixSparse::setRowToConstant(int irow, double value)
{
  int ncols = getNCols();
  if (getFlagMatrixCheck())
  {
    if (!_isRowValid(irow)) return;
  }
  for (int icol = 0; icol < ncols; icol++)
    eigenMat().coeffRef(irow, icol) = value;
}

void MatrixSparse::setDiagonal(const VectorDouble& tab)
{
  if (!isSquare())
    my_throw("This function is only valid for Square matrices");
  if (getFlagMatrixCheck())
  {
    if (!_isRowSizeConsistent(tab)) return;
  }

  Eigen::Map<const Eigen::VectorXd> vecm(tab.data(), tab.size());
  eigenMat() = vecm.asDiagonal();
}

void MatrixSparse::setDiagonalToConstant(double value)
{
  if (!isSquare())
    my_throw("This function is only valid for Square matrices");

  VectorDouble vec(getNRows(), value);
  Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
  eigenMat() = vecm.asDiagonal();
}

/*! Gets the value for rank 'rank' */
double MatrixSparse::_getValueByRank(int rank) const
{
  DECLARE_UNUSED(rank);
  _forbiddenForSparse("_getValueByRank");
  return TEST;
}

double& MatrixSparse::_getValueRef(int irow, int icol)
{
  DECLARE_UNUSED(irow);
  DECLARE_UNUSED(icol);
  _forbiddenForSparse("_getValueRef");
  return AMatrix::_getValueRef(irow, icol);
}

void MatrixSparse::_setValueByRank(int rank, double value)
{
  DECLARE_UNUSED(rank);
  DECLARE_UNUSED(value);
  _forbiddenForSparse("_setValueByRank");
}

void MatrixSparse::setValue(int irow, int icol, double value)
{
  if (getFlagMatrixCheck() && !_isIndexValid(irow, icol)) return;
  eigenMat().coeffRef(irow, icol) = value;
}

void MatrixSparse::updValue(int irow,
                            int icol,
                            const EOperator& oper,
                            double value)
{
  if (getFlagMatrixCheck() && !_isIndexValid(irow, icol)) return;
  double newval                   = modifyOperator(oper, eigenMat().coeff(irow, icol), value);
  eigenMat().coeffRef(irow, icol) = newval;
}

int MatrixSparse::_getMatrixPhysicalSize() const
{
  return eigenMat().nonZeros();
}

/**
 * @param value Constant value used for filling 'this'
 */
void MatrixSparse::fill(double value)
{
  int nrow = getNRows();
  int ncol = getNCols();
  NF_Triplet NF_T;
  for (int irow = 0; irow < nrow; irow++)
    for (int icol = 0; icol < ncol; icol++)
      NF_T.add(irow, icol, value);

  resetFromTriplet(NF_T);
}

/*! Multiply a Matrix row-wise */
void MatrixSparse::multiplyRow(const VectorDouble& vec)
{
  for (int k = 0; k < eigenMat().outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), k); it; ++it)
      it.valueRef() *= vec[it.row()];
}

/*! Multiply a Matrix column-wise */
void MatrixSparse::multiplyColumn(const VectorDouble& vec)
{
  for (int k = 0; k < eigenMat().outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), k); it; ++it)
      it.valueRef() *= vec[it.col()];
}

/*! Divide a Matrix row-wise */
void MatrixSparse::divideRow(const VectorDouble& vec)
{
  for (int k = 0; k < eigenMat().outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), k); it; ++it)
      it.valueRef() /= vec[it.row()];
}

/*! Divide a Matrix column-wise */
void MatrixSparse::divideColumn(const VectorDouble& vec)
{
  for (int k = 0; k < eigenMat().outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), k); it; ++it)
      it.valueRef() /= vec[it.col()];
}

void MatrixSparse::prodMatVecInPlaceC(constvect x, vect res, bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x.data(), x.size());
  Eigen::Map<Eigen::VectorXd> ym(res.data(), res.size());
  if (transpose)
    ym = eigenMat().transpose() * xm;
  else
    ym = eigenMat() * xm;
}

/*! Perform y += 'this' %*% x */
void MatrixSparse::addProdMatVecInPlaceToDest(const constvect in,
                                              vect out,
                                              bool transpose) const
{
  Eigen::Map<const Eigen::VectorXd> inm(in.data(), in.size());
  Eigen::Map<Eigen::VectorXd> outm(out.data(), out.size());
  if (transpose)
    outm += eigenMat().transpose() * inm;
  else
    outm += eigenMat() * inm;
}

/**
 * Filling the matrix with an array of values
 * Note that this array is ALWAYS dimensioned to the total number
 * of elements in the matrix.
 * Kept for compatibility with old code where matrix contents was stored as
 * a double* array
 * @param values Input array (Dimension: nrow * ncol)
 * @param byCol true for Column Major; false for Row Major
 */
#ifndef SWIG
void MatrixSparse::_setValues(const double* values, bool byCol)
{
  if (byCol)
  {
    Eigen::Map<const Eigen::MatrixXd> temp(values, getNRows(), getNCols());
    eigenMat() = temp.sparseView(1., EPSILON10);
  }
  else
  {
    Eigen::Map<const Eigen::MatrixXd> temp(values, getNCols(), getNRows());
    eigenMat() = temp.transpose().sparseView(1., EPSILON10);
  }
}
#endif

MatrixSparse* MatrixSparse::create(const MatrixSparse* mat)
{
  return new MatrixSparse(*mat);
}

MatrixSparse* MatrixSparse::create(int nrow, int ncol)
{
  return new MatrixSparse(nrow, ncol);
}

MatrixSparse* MatrixSparse::createFromTriplet(const NF_Triplet& NF_T,
                                              int nrow,
                                              int ncol,
                                              int nrowmax)
{
  // If 'nrow' a  nd 'ncol' are not defined, derive them from NF_T
  if (nrow <= 0 || ncol <= 0)
  {
    nrow = NF_T.getNRows() + 1;
    ncol = NF_T.getNCols() + 1;
  }
  MatrixSparse* mat = new MatrixSparse(nrow, ncol, nrowmax);
  mat->resetFromTriplet(NF_T);

  return mat;
}

MatrixSparse* MatrixSparse::Identity(int nrow, double value)
{
  MatrixSparse* mat = new MatrixSparse(nrow, nrow);
  for (int i = 0; i < nrow; i++)
  {
    mat->eigenMat().coeffRef(i, i) += value;
  }
  return mat;
}

MatrixSparse* MatrixSparse::addMatMat(const MatrixSparse* x,
                                      const MatrixSparse* y,
                                      double cx,
                                      double cy)
{
  MatrixSparse* mat = new MatrixSparse(x->getNRows(), x->getNCols(), -1);
  mat->eigenMat()   = cx * x->eigenMat() + cy * y->eigenMat();
  return mat;
}

MatrixSparse* MatrixSparse::diagVec(const VectorDouble& vec)
{
  int size          = (int)vec.size();
  MatrixSparse* mat = new MatrixSparse(size, size);

  mat->setDiagonal(vec);
  return mat;
}

MatrixSparse* MatrixSparse::diagConstant(int number, double value)
{
  MatrixSparse* mat = new MatrixSparse(number, number);
  mat->setDiagonalToConstant(value);
  return mat;
}

/**
 * Construct a sparse matrix with the diagonal of 'A', where each element is transformed
 * @param A    Input sparse matrix
 * @param oper_choice: Operation on the diagonal term (see Utilities::operate_XXX)
 * @return
 */
MatrixSparse* MatrixSparse::diagMat(MatrixSparse* A, int oper_choice)
{
  if (!A->isSquare())
  {
    messerr("This method requires the matrix 'A' to be square");
    return nullptr;
  }

  VectorDouble diag = A->getDiagonal();
  VectorHelper::transformVD(diag, oper_choice);
  return MatrixSparse::diagVec(diag);
  return MatrixSparse::diagVec(diag);
}

bool MatrixSparse::_isElementPresent(int irow, int icol) const
{
  for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), icol); it; ++it)
  {
    if (it.row() == irow) return true;
  }
  return false;
}

void MatrixSparse::addValue(int row, int col, double value)
{
  if (ABS(value) <= EPSILON10) return;
  eigenMat().coeffRef(row, col) += value;
}

double MatrixSparse::getValue(int row, int col) const
{
  if (getFlagMatrixCheck() && !_isIndexValid(row, col)) return TEST;
  return eigenMat().coeff(row, col);
}

double MatrixSparse::L1Norm() const
{
  return (Eigen::RowVectorXd::Ones(eigenMat().rows()) * eigenMat().cwiseAbs()).maxCoeff();
}

void MatrixSparse::getStats(int* nrows, int* ncols, int* count, double* percent) const
{
  *nrows   = getNRows();
  *ncols   = getNCols();
  *count   = eigenMat().nonZeros();
  *percent = 0.;
  if ((*nrows) > 0 && (*ncols) > 0)
    (*percent) = ((100. * (double)(*count)) / ((double)(*nrows) * (double)(*ncols)));
}

VectorDouble MatrixSparse::extractDiag(int oper_choice) const
{
  VectorDouble diag(std::min(getNCols(), getNRows()));
  Eigen::Map<Eigen::VectorXd> ym(diag.data(), diag.size());
  ym = eigenMat().diagonal();
  VH::transformVD(diag, oper_choice);
  return diag;
}

int MatrixSparse::addVecInPlaceEigen(const Eigen::Map<const Eigen::VectorXd>& xm,
                                     Eigen::Map<Eigen::VectorXd>& ym) const
{
  ym = eigenMat() * xm + ym;
  return 0;
}

int MatrixSparse::addVecInPlace(const constvect xm, vect ym) const
{
  Eigen::Map<const Eigen::VectorXd> xmm(xm.data(), xm.size());
  Eigen::Map<Eigen::VectorXd> ymm(ym.data(), ym.size());
  ymm = eigenMat() * xmm + ymm;
  return 0;
}

int MatrixSparse::addVecInPlaceVD(const VectorDouble& x, VectorDouble& y) const
{
  Eigen::Map<const Eigen::VectorXd> xm(x.data(), x.size());
  Eigen::Map<Eigen::VectorXd> ym(y.data(), y.size());
  ym = eigenMat() * xm + ym;
  return 0;
}

void MatrixSparse::setConstant(double value)
{
  for (int k = 0; k < eigenMat().outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), k); it; ++it)
      it.valueRef() = value;
}

int MatrixSparse::scaleByDiag()
{
  VectorDouble diag = extractDiag(-1);
  Eigen::Map<const Eigen::VectorXd> ym(diag.data(), diag.size());
  eigenMat() = ym.asDiagonal() * eigenMat();
  return 0;
}

/**
 *
 * @param v Add a scalar value to all terms of the current matrix
 */
void MatrixSparse::addScalar(double v)
{
  if (isZero(v)) return;
  for (int k = 0; k < eigenMat().outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), k); it; ++it)
      it.valueRef() += v;
}

/**
 *
 * @param v Add constant value to the diagonal of the current Matrix
 */
void MatrixSparse::addScalarDiag(double v)
{
  if (isZero(v)) return;

  for (int k = 0; k < eigenMat().outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), k); it; ++it)
    {
      if (it.col() == it.row())
        it.valueRef() += v;
    }
}

/**
 *
 * @param v Multiply all the terms of the matrix by the scalar 'v'
 */
void MatrixSparse::prodScalar(double v)
{
  if (isOne(v)) return;
  for (int k = 0; k < eigenMat().outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), k); it; ++it)
      it.valueRef() *= v;
}

void MatrixSparse::_addProdMatVecInPlacePtr(const double* x, double* y, bool transpose) const
{
  if (transpose)
  {
    Eigen::Map<const Eigen::VectorXd> xm(x, getNRows());
    Eigen::Map<Eigen::VectorXd> ym(y, getNCols());
    ym += eigenMat().transpose() * xm;
  }
  else
  {
    Eigen::Map<const Eigen::VectorXd> xm(x, getNCols());
    Eigen::Map<Eigen::VectorXd> ym(y, getNRows());
    ym += eigenMat() * xm;
  }
}

void MatrixSparse::_addProdVecMatInPlacePtr(const double* x, double* y, bool transpose) const
{
  if (transpose)
  {
    Eigen::Map<const Eigen::VectorXd> xm(x, getNCols());
    Eigen::Map<Eigen::VectorXd> ym(y, getNRows());
    ym += xm.transpose() * eigenMat().transpose();
  }
  else
  {
    Eigen::Map<const Eigen::VectorXd> xm(x, getNRows());
    Eigen::Map<Eigen::VectorXd> ym(y, getNCols());
    ym += xm.transpose() * eigenMat();
  }
}

/**
 * Store the product of 'x' by 'y' in this
 * @param x First Matrix
 * @param y Second matrix
 * @param transposeX True if First matrix is transposed
 * @param transposeY True if Second matrix is transposed
 */
void MatrixSparse::prodMatMatInPlace(const AMatrix* x,
                                     const AMatrix* y,
                                     bool transposeX,
                                     bool transposeY)
{
  if (getFlagMatrixCheck() &&
      !_isMatrixCompatible("MatrixSparse::prodMatMatInPlace",
                           x, 0, transposeX,
                           y, 0, transposeY)) return;

  const MatrixSparse* xm = dynamic_cast<const MatrixSparse*>(x);
  const MatrixSparse* ym = dynamic_cast<const MatrixSparse*>(y);
  if (xm == nullptr || ym == nullptr)
  {
    AMatrix::prodMatMatInPlace(x, y, transposeX, transposeY);
  }
  else
  {
    if (transposeX)
    {
      if (transposeY)
        eigenMat() = xm->eigenMat().transpose() * ym->eigenMat().transpose();
      else
        eigenMat() = xm->eigenMat().transpose() * ym->eigenMat();
    }
    else
    {
      if (transposeY)
        eigenMat() = xm->eigenMat() * ym->eigenMat().transpose();
      else
        eigenMat() = xm->eigenMat() * ym->eigenMat();
    }
  }
}

MatrixSparse* prodNormMatMat(const MatrixSparse* a,
                             const MatrixSparse* m,
                             bool transpose)
{
  int nrow          = (transpose) ? a->getNCols() : a->getNRows();
  int ncol          = (transpose) ? a->getNRows() : a->getNCols();
  MatrixSparse* mat = new MatrixSparse(nrow, ncol);
  mat->prodNormMatMatInPlace(a, m, transpose);
  return mat;
}

MatrixSparse* prodNormMatVec(const MatrixSparse* a, const VectorDouble& vec, bool transpose)
{
  int nsym          = (transpose) ? a->getNCols() : a->getNRows();
  MatrixSparse* mat = new MatrixSparse(nsym, nsym);
  mat->prodNormMatVecInPlace(a, vec, transpose);
  return mat;
}

MatrixSparse* prodNormMat(const MatrixSparse* a, bool transpose)
{
  int nsym          = (transpose) ? a->getNCols() : a->getNRows();
  MatrixSparse* mat = new MatrixSparse(nsym, nsym);
  mat->prodNormMatInPlace(a, transpose);
  return mat;
}

MatrixSparse* prodNormDiagVec(const MatrixSparse* a,
                              const VectorDouble& vec,
                              int oper_choice)
{
  int nrow          = a->getNRows();
  int ncol          = a->getNCols();
  MatrixSparse* mat = new MatrixSparse(nrow, ncol, -1);

  // Perform the transformation of the input vector
  VectorDouble vecp = vec;
  VH::transformVD(vecp, oper_choice);

  Eigen::Map<const Eigen::VectorXd> vecm(vecp.data(), vecp.size());
  auto diag       = vecm.asDiagonal();
  mat->eigenMat() = (diag * a->eigenMat() * diag);
  return mat;
}

/**
 * Perform: 'this' = diag('vec') %*% 'A' %*% diag('vec')
 * @param vec  Input Vector
 * @param oper_choice Type of transformation
 */
void MatrixSparse::prodNormDiagVecInPlace(const VectorDouble& vec, int oper_choice)
{
  if (!isSquare())
  {
    messerr("This method is limited to square matrices");
    return;
  }
  if (getNRows() != (int)vec.size())
  {
    messerr("Matrix dimension (%d) does not match vector dimension (%d)",
            getNRows(), (int)vec.size());
    return;
  }

  // Perform the transformation of the input vector
  VectorDouble vecp = vec;
  VH::transformVD(vecp, oper_choice);

  Eigen::Map<const Eigen::VectorXd> vecm(vecp.data(), vecp.size());
  auto diag  = vecm.asDiagonal();
  eigenMat() = diag * eigenMat() * diag;
}

void MatrixSparse::prodNormMatVecInPlace(const AMatrix* a,
                                         const VectorDouble& vec,
                                         bool transpose)
{
  if (getFlagMatrixCheck())
  {
    // Note that 'vec' is not tested for compatibility:
    // it is a vector but used as a diagonal of a square matrix
    if (!_isMatrixCompatible("MatrixSparse::prodNormMatVecInPlace",
                             a, 0, transpose,
                             a, 0, !transpose)) return;
  }

  const MatrixSparse* am = dynamic_cast<const MatrixSparse*>(a);
  if (am == nullptr)
  {
    AMatrix::prodNormMatVecInPlace(a, vec, transpose);
  }
  else
  {
    if (transpose)
    {
      Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
      eigenMat() = am->eigenMat().transpose() * vecm.asDiagonal() * am->eigenMat();
    }
    else
    {
      Eigen::Map<const Eigen::VectorXd> vecm(vec.data(), vec.size());
      eigenMat() = am->eigenMat() * vecm.asDiagonal() * am->eigenMat().transpose();
    }
  }
}

void MatrixSparse::prodNormMatInPlace(const AMatrix* a, bool transpose)
{
  if (getFlagMatrixCheck() &&
      !_isMatrixCompatible("MatrixSparse::prodNormMatInPlace",
                           a, 0, transpose,
                           a, 0, !transpose)) return;

  const MatrixSparse* am = dynamic_cast<const MatrixSparse*>(a);
  if (am == nullptr)
  {
    AMatrix::prodNormMatInPlace(a, transpose);
  }
  else
  {
    if (transpose)
    {
      eigenMat() = am->eigenMat().transpose() * am->eigenMat();
    }
    else
    {
      eigenMat() = am->eigenMat() * am->eigenMat().transpose();
    }
  }
}

void MatrixSparse::prodNormMatMatInPlace(const AMatrix* a,
                                         const AMatrix* m,
                                         bool transpose)
{
  if (getFlagMatrixCheck() &&
      !_isMatrixCompatible("MatrixSparse::prodNormMatMatInPlace",
                           a, 0, transpose,
                           m, 0, false,
                           a, 0, !transpose)) return;

  const MatrixSparse* am = dynamic_cast<const MatrixSparse*>(a);
  const MatrixSparse* mm = dynamic_cast<const MatrixSparse*>(m);
  if (am == nullptr || mm == nullptr)
  {
    AMatrix::prodNormMatMatInPlace(a, m, transpose);
  }
  else
  {
    if (transpose)
    {
      eigenMat() = (am->eigenMat().transpose() * mm->eigenMat()) * am->eigenMat();
    }
    else
    {
      eigenMat() = (am->eigenMat() * mm->eigenMat()) * am->eigenMat().transpose();
    }
  }
}

/*!
 * Updates the current Matrix as a linear combination of matrices as follows:
 *  this <- cx * this + cy * y
 * @param cx Coefficient applied to the current Matrix
 * @param cy Coefficient applied to the Matrix  'y'
 * @param y Second Matrix in the Linear combination
 */
void MatrixSparse::addMatInPlace(const MatrixSparse& y, double cx, double cy)
{
  if (getFlagMatrixCheck() &&
      !isSameSize(y)) return;

  eigenMat() = cx * eigenMat() + cy * y.eigenMat();
}

int MatrixSparse::_invert()
{
  if (!isSquare())
    my_throw("Invert method is restricted to Square matrices");
  int n = getNCols();
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(eigenMat());
  Eigen::SparseMatrix<double> I(n, n);
  I.setIdentity();
  eigenMat() = solver.solve(I);

  return 0;
}

int MatrixSparse::_solve(const VectorDouble& b, VectorDouble& x) const
{
  if (!isSquare())
    my_throw("Invert method is limited to Square Matrices");
  if ((int)b.size() != getNRows() || (int)x.size() != getNRows())
    my_throw("b' and 'x' should have the same dimension as the Matrix");

  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  Eigen::Map<const Eigen::VectorXd> bm(b.data(), getNCols());
  Eigen::Map<Eigen::VectorXd> xm(x.data(), getNRows());
  xm = solver.compute(eigenMat()).solve(bm);

  return 0;
}

String MatrixSparse::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  sstr << AMatrix::toString(strfmt) << std::endl;
  return sstr.str();
}

void MatrixSparse::_allocate(int nrow, int ncol, int ncolmax)
{
  eigenMat() = Eigen::SparseMatrix<double>(nrow, ncol);
  _setNCols(ncol);
  _setNRows(nrow);

  if (ncolmax > 0)
  {
    eigenMat().reserve(Eigen::VectorXi::Constant(nrow, ncolmax));
  }
  _nColMax = ncolmax;
}

/**
 * This strange function instantiate a sparse matrix with given dimensions
 * filled with zeroes. It should be an empty matrix... But this does not make sense.
 * Therefore it is created by setting a single element at the lower bottom size of
 * the matrix ... filled with a zero.
 */
void MatrixSparse::_allocate()
{
  eigenMat() = Eigen::SparseMatrix<double>(getNRows(), getNCols());
  {
    eigenMat() = Eigen::SparseMatrix<double>(getNRows(), getNCols());

    if (_nColMax > 0)
    {
      eigenMat().reserve(Eigen::VectorXi::Constant(getNCols(), _nColMax));
    }
    if (isMultiThread()) omp_set_num_threads(getMultiThread());
  }
}

void MatrixSparse::_deallocate()
{
  eigenMat().data().squeeze();
  {
    eigenMat().data().squeeze();
  }
}

void MatrixSparse::_forbiddenForSparse(const String& func)
{
  messerr("Problem with Function: %s", func.c_str());
  messerr("This function is not available in Sparse Matrix");
}

void MatrixSparse::dumpElements(const String& title, int ifrom, int ito)
{
  DECLARE_UNUSED(title);
  DECLARE_UNUSED(ifrom);
  DECLARE_UNUSED(ito);
  messerr("This method is not implemented for Sparse Matrix");
}

/**
 * From a matrix of any type, creates the triplet
 * (specific format for creating efficiently a Sparse matrix)
 * It only takes the only non-zero elements of the matrix
 */
NF_Triplet MatrixSparse::getMatrixToTriplet(int shiftRow, int shiftCol) const
{
  return NF_Triplet::createFromEigen(eigenMat(), shiftRow, shiftCol);
}

void MatrixSparse::_clear()
{
  _setNRows(0);
  _setNCols(0);
  _allocate();
}

int MatrixSparse::_getIndexToRank(int irow, int icol) const
{
  DECLARE_UNUSED(irow);
  DECLARE_UNUSED(icol);
  _forbiddenForSparse("_getIndexToRank");
  return ITEST;
}

MatrixSparse* createFromAnyMatrix(const AMatrix* matin)
{
  return MatrixSparse::createFromTriplet(matin->getMatrixToTriplet(),
                                         matin->getNRows(),
                                         matin->getNCols(),
                                         -1);
}

int MatrixSparse::_eigen_findColor(int imesh,
                                   int ncolor,
                                   VectorInt& colors,
                                   VectorInt& temp) const
{
  temp.fill(0);

  /* Checks the colors of the connected nodes */

  for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), imesh); it; ++it)
  {
    if (isZero(it.value())) continue;
    int irow = it.row();
    if (!IFFFF(colors[irow])) temp[colors[irow] - 1]++;
  }

  /* Look for a free color */

  for (int j = 0; j < ncolor; j++)
  {
    if (temp[j] == 0) return (j + 1);
  }
  return (-1);
}

VectorInt MatrixSparse::colorCoding() const
{
  int next_col = 0;
  int ncol     = 0;
  int nmesh    = getNCols();

  /* Core allocation */

  VectorInt colors(nmesh, ITEST);
  VectorInt temp(nmesh);

  /* Loop on the nodes of the mesh */

  for (int imesh = 0; imesh < nmesh; imesh++)
  {
    next_col = _eigen_findColor(imesh, ncol, colors, temp);
    next_col = _eigen_findColor(imesh, ncol, colors, temp);

    if (next_col < 0)
    {
      ncol++;
      colors[imesh] = ncol;
    }
    else
    {
      colors[imesh] = next_col;
    }
  }
  return colors;
}

void MatrixSparse::glueInPlace(MatrixSparse* A1,
                               const MatrixSparse* A2,
                               bool flagShiftRow,
                               bool flagShiftCol)
{
  int shiftRow  = (flagShiftRow) ? A1->getNRows() : 0;
  int shiftCol  = (flagShiftCol) ? A1->getNCols() : 0;
  NF_Triplet T1 = A1->getMatrixToTriplet();
  NF_Triplet T2 = A2->getMatrixToTriplet(shiftRow, shiftCol);

  // Concatenate the two triplet lists
  T1.appendInPlace(T2);

  // Create the new matrix from the resulting triplet list
  int nrow = (flagShiftRow) ? A1->getNRows() + A2->getNRows() : MAX(A1->getNRows(), A2->getNRows());
  int ncol = (flagShiftCol) ? A1->getNCols() + A2->getNCols() : MAX(A1->getNCols(), A2->getNCols());

  A1->resetFromTriplet(T1);
  A1->_setNRows(nrow);
  A1->_setNCols(ncol);
}

MatrixSparse* MatrixSparse::glue(const MatrixSparse* A1,
                                 const MatrixSparse* A2,
                                 bool flagShiftRow,
                                 bool flagShiftCol)
{
  int shiftRow = (flagShiftRow) ? A1->getNRows() : 0;
  int shiftCol = (flagShiftCol) ? A1->getNCols() : 0;

  // Create the two triplet lists
  NF_Triplet T1 = A1->getMatrixToTriplet();
  NF_Triplet T2 = A2->getMatrixToTriplet(shiftRow, shiftCol);

  // Concatenate the two triplet lists
  T1.appendInPlace(T2);

  // Create the new matrix from the resulting triplet list
  int nrow = (flagShiftRow) ? A1->getNRows() + A2->getNRows() : MAX(A1->getNRows(), A2->getNRows());
  int ncol = (flagShiftCol) ? A1->getNCols() + A2->getNCols() : MAX(A1->getNCols(), A2->getNCols());

  return MatrixSparse::createFromTriplet(T1, nrow, ncol, -1);
  return MatrixSparse::createFromTriplet(T1, nrow, ncol, -1);
}

/* Extract a sparse sub-matrix */
/* 'rank_rows' and 'rank_cols' must have same dimension as C */
/* The arrays 'rank_rows' and 'rank_cols' may be absent */
/* Their value gives the rank of the saved element or -1 */
MatrixSparse* MatrixSparse::extractSubmatrixByRanks(const VectorInt& rank_rows,
                                                    const VectorInt& rank_cols) const
{
  int old_row, old_col, new_row, new_col;

  NF_Triplet NF_Tin = getMatrixToTriplet();
  NF_Triplet NF_Tout;

  /* Fill the new sparse triplet */

  for (int i = 0; i < NF_Tin.getNElements(); i++)
  {
    old_row = NF_Tin.getRow(i);
    old_col = NF_Tin.getCol(i);
    new_row = (!rank_rows.empty()) ? rank_rows[old_row] : old_row;
    new_col = (!rank_cols.empty()) ? rank_cols[old_col] : old_col;
    if (new_row < 0 || new_col < 0) continue;
    NF_Tout.add(new_row, new_col, NF_Tin.getValue(i));
  }

  return MatrixSparse::createFromTriplet(NF_Tout);
}

/* Extract a sparse submatrix */
/* The array 'colors' has the same dimension as C */
/* The element of 'C' must be kept if: */
/* - the color of its row number is equal to 'ref_color' if 'row_ok'==TRUE */
/*   or different if 'row_ok'== FALSE */
/* and if */
/* - the color of its column number is equal to 'ref-color' if 'col_ok'==TRUE*/
/*   or different if 'col_ok'== FALSE */
MatrixSparse* MatrixSparse::extractSubmatrixByColor(const VectorInt& colors,
                                                    int ref_color,
                                                    bool row_ok,
                                                    bool col_ok)
{
  /* Convert the contents of the sparse matrix into columns */

  NF_Triplet NF_Tin = getMatrixToTriplet();

  /* Initialize the output matrix */

  NF_Triplet NF_Tout;

  /* Core allocation */

  int n = getNCols();
  VectorInt u_row(n);
  VectorInt u_col(n);

  int ir = 0;
  for (int i = 0; i < n; i++)
  {
    u_row[i] = -1;
    if (row_ok && colors[i] != ref_color) continue;
    if (!row_ok && colors[i] == ref_color) continue;
    u_row[i] = ir++;
  }

  int ic = 0;
  for (int i = 0; i < n; i++)
  {
    u_col[i] = -1;
    if (col_ok && colors[i] != ref_color) continue;
    if (!col_ok && colors[i] == ref_color) continue;
    u_col[i] = ic++;
  }

  /* Fill the new sparse triplet */

  for (int i = 0; i < NF_Tin.getNElements(); i++)
  {
    ir = u_row[NF_Tin.getRow(i)];
    ic = u_col[NF_Tin.getCol(i)];
    if (ir < 0 || ic < 0) continue;
    NF_Tout.add(ir, ic, NF_Tin.getValue(i));
  }

  return MatrixSparse::createFromTriplet(NF_Tout, 0, 0, -1);
}

void MatrixSparse::gibbs(int iech,
                         const VectorDouble& zcur,
                         double* yk,
                         double* sk)
{
  *yk = 0.;
  for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), iech); it;
       ++it)
  {
    double coeff = it.valueRef();
    if (ABS(coeff) <= 0.) continue;
    int jech = it.row();
    {
      if (iech == jech)
        *sk = coeff;
      else
        *yk -= coeff * zcur[jech];
    }
  }

  // Returned arguments
  (*yk) /= (*sk);
  (*sk) = sqrt(1. / (*sk));
}

int MatrixSparse::_addToDest(const constvect inv, vect outv) const
{
  Eigen::Map<const Eigen::VectorXd> inm(inv.data(), inv.size());
  Eigen::Map<Eigen::VectorXd> outm(outv.data(), outv.size());
  outm += eigenMat() * inm;
  return 0;
}

void MatrixSparse::setDiagonal(const constvect tab)
{
  Eigen::Map<const Eigen::VectorXd> tabm(tab.data(), tab.size());
  setDiagonal(tabm);
}

void MatrixSparse::setDiagonal(const Eigen::Map<const Eigen::VectorXd>& tab)
{
  eigenMat() = tab.asDiagonal();
}

/*! Extract a Row */
MatrixSparse* MatrixSparse::getRowAsMatrixSparse(int irow, double coeff) const
{
  int ncols         = getNCols();
  MatrixSparse* res = new MatrixSparse(1, ncols);

  // The input sparse matrix being symmetrical, we benefit from its
  // column-major storage (setting icol = irow)
  int icol = irow;
  for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), icol); it; ++it)
    res->eigenMat().coeffRef(0, it.row()) = coeff * it.value();

  return res;
}

/*! Extract a Column */
MatrixSparse* MatrixSparse::getColumnAsMatrixSparse(int icol, double coeff) const
{
  int nrows         = getNRows();
  MatrixSparse* res = new MatrixSparse(nrows, 1);

  for (Eigen::SparseMatrix<double>::InnerIterator it(eigenMat(), icol); it; ++it)
    res->eigenMat().coeffRef(it.row(), 0) = coeff * it.value();

  return res;
}

///////////////Not exported //////////

Eigen::SparseMatrix<double> AtMA(const Eigen::SparseMatrix<double>& A,
                                 const Eigen::SparseMatrix<double>& M)
{
  return A.transpose() * M * A;
}

int MatrixSparse::forwardLU(const VectorDouble& b, VectorDouble& x, bool flagLower) const
{
  Eigen::Map<const Eigen::VectorXd> bm(b.data(), b.size());
  Eigen::Map<Eigen::VectorXd> xm(x.data(), x.size());

  if (!flagLower)
  {
    const Eigen::SparseMatrix<double>& Lx = eigenMat().transpose();
    xm                                    = Lx.triangularView<Eigen::Upper>().solve(bm);
  }
  else
  {
    const Eigen::SparseMatrix<double>& Lx = eigenMat();
    xm                                    = Lx.triangularView<Eigen::Lower>().solve(bm);
  }
  return 0;
}

void MatrixSparse::forceDimension(int maxRows, int maxCols)
{
  // Redimensionner les matrices si nécessaire
  if (eigenMat().rows() < maxRows || eigenMat().cols() < maxCols)
  {
    eigenMat().conservativeResize(maxRows, maxCols);
    eigenMat().insert(maxRows - 1, maxCols - 1) = 0.0; // Élément fictif
  }
}
} // namespace gstlrn