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
#include "Estimation/Likelihood.hpp"
#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Db/RankHandler.hpp"
#include "Estimation/ALikelihood.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/ModelGeneric.hpp"
#include "Space/SpacePoint.hpp"
#include "Stats/Classical.hpp"
#include "Tree/Ball.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
Likelihood::Likelihood(ModelGeneric* model,
                       const Db* db,
                       bool reml)
  : ALikelihood(model, db, reml)
{
  setAuthorizedAnalyticalGradients(true);
}

Likelihood::Likelihood(const Likelihood& r)
  : ALikelihood(r)
{
}

Likelihood& Likelihood::operator=(const Likelihood& r)
{
  if (this != &r)
  {
    ALikelihood::operator=(r);
  }
  return *this;
}

Likelihood::~Likelihood()
{
}

double logLikelihood(const Db* db,
                     ModelGeneric* model,
                     bool verbose)
{
  Likelihood* vec = Likelihood::createForOptim(model, db);
  double result   = vec->computeLogLikelihood(verbose);
  delete vec;
  return result;
}

Likelihood* Likelihood::createForOptim(ModelGeneric* model,
                                       const Db* db,
                                       bool reml)
{
  auto* vec            = new Likelihood(model, db, reml);
  MatrixSymmetric vars = dbVarianceMatrix(db);
  double hmax          = db->getExtensionDiagonal();
  vec->setEnvironment(vars, hmax);
  vec->init();
  return vec;
}

void Likelihood::_computeCm1X()
{
  if (_covChol.solveMatrix(_X, _Cm1X))
  {
    messerr("Problem when solving a Linear System after Cholesky decomposition");
  }
}

void Likelihood::_computeCm1Y()
{
  _Cm1Y.resize(_Y.size());
  if (_covChol.solve(_Y, _Cm1Y))
  {
    messerr("Error when calculating Cm1Z");
  }
}

double Likelihood::_computeLogDet() const
{
  return _covChol.computeLogDeterminant();
}

void Likelihood::_updateModel(bool verbose)
{
  DECLARE_UNUSED(verbose);
  _model->evalCovMatSymInPlace(_cov, _db);
  _covChol.setMatrix(&_cov);
}

void Likelihood::evalGrad(vect res)
{

  _temp.resize(_Y.size());
  _gradCovMatTimesInvCov.resize(_Y.size(), _Y.size());
  auto invcov = _covChol.inverse();
  RankHandler rkh(_db);
  rkh.defineSampleRanks();
  auto gradcov = _model->getGradients();
  _gradCovMat.resize(_Y.size(), _Y.size());
  CholeskyDense XtCm1XChol;
  MatrixSymmetric invXtCm1X;
  if (_reml && _model->getNDriftEquation() > 0)
  {
    XtCm1XChol.setMatrix(&_XtCm1X);
    invXtCm1X = XtCm1XChol.inverse();
  }
  for (size_t iparam = 0; iparam < gradcov.size(); iparam++)
  {
    _fillGradCovMat(rkh, gradcov[iparam]);
    _gradCovMat.prodMatVecInPlace(_Cm1Y, _temp);
    double dquad = -VH::innerProduct(_Cm1Y, _temp);
    res[iparam] = 0.0;
    if (_reml && _model->getNDriftEquation() > 0)
    {
      MatrixSymmetric temp(_X.getNCols());
      temp.prodNormMatMatInPlace(&_Cm1X, &_gradCovMat, true);
      double dlogdetreml = MatrixDense::traceProd(invXtCm1X, temp);
      res[iparam] += 0.5 * dlogdetreml;
    }
    double dlogdet = MatrixDense::traceProd(invcov, _gradCovMat); // Warning: _gradCovMat is modified so the line
    // has to be after _gradCovMat.prodMatVecInPlace(_Cm1Y, _temp);
    res[iparam] += 0.5 * (dlogdet + dquad);

  }
}

void Likelihood::_fillGradCovMat(RankHandler& rkh, covmaptype& gradcov)
{
  int icur, jcur = 0;

  SpacePoint p1, p2;
  rkh.defineSampleRanks();

  for (size_t jvar = 0; (int)jvar < _model->getNVar(); jvar++)
  {
    auto indsj = rkh.getSampleRanksByVariable(jvar);

    for (auto& j: indsj)
    {
      icur = 0;
      _db->getSampleAsSPInPlace(p1, j);

      for (size_t ivar = 0; (int)ivar < _model->getNVar(); ivar++)
      {
        auto indsi = rkh.getSampleRanksByVariable(ivar);
        for (auto& i: indsi)
        {
          _db->getSampleAsSPInPlace(p2, i);

          _gradCovMat.setValue(icur, jcur, gradcov(p1, p2, ivar, jvar, nullptr));
          icur++;
        }
      }
      jcur++;
    }
  }
}
}