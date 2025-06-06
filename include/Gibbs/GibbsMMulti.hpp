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

#include "Gibbs/GibbsMulti.hpp"
#include "LinearOp/CholeskySparse.hpp"

class MatrixSparse;
class Db;
class Model;

class GSTLEARN_EXPORT GibbsMMulti: public GibbsMulti
{
public:
  GibbsMMulti();
  GibbsMMulti(Db* db, Model* model);
  GibbsMMulti(const GibbsMMulti &r);
  GibbsMMulti& operator=(const GibbsMMulti &r);
  virtual ~GibbsMMulti();

  void update(VectorVectorDouble &y, int isimu, int ipgs, int iter) override;
  int covmatAlloc(bool verbose, bool verboseTimer = false) override;

  void setEps(double eps) { _eps = eps; }
  void cleanup() override;

  void setFlagStoreInternal(bool flagStoreInternal) ;

private:
  int  _getNVar() const;
  int  _getSize() const;
  void _storeWeights(int icol);
  void _storeWeightsMS(int icol, NF_Triplet& NF_T);
  void _getWeights(int icol) const;
  void _calculateWeights(int icol,
                         VectorDouble &b,
                         double tol = EPSILON3) const;
  int  _storeAllWeights(bool verbose = false);
  int  _getSizeOfWeights(const VectorDouble& weights) const;
  double _getEstimate(int ipgs, int icol, const VectorVectorDouble &y) const;
  void _allocate();
  double _getVariance(int icol) const;
  int _getColumn(int iact, int ivar) const;
  void _splitCol(int icol, int *iact, int *ivar) const;
  void _updateStatWeights(int* nzero);

private:
  MatrixSparse*   _Cmat;
  CholeskySparse* _CmatChol;
  double          _eps;
  bool            _flagStoreInternal;

  VectorVectorDouble _areas;
  MatrixSparse*      _matWgt;

  mutable VectorDouble _weights;
};
