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
#include "Covariances/CovPotential.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

#define TR(i, j) (_Tr[3 * (i) + (j)])

namespace gstlrn
{
CovPotential::CovPotential(const CovAniso& cova)
  : ACov()
  , _nVar(0)
  , _covRef(cova)
  , _launchCalculations(true)
  , _flagGradient(false)
  , _covpp(0.)
  , _covGp()
  , _covGG()
{
  setContext(cova.getContext());
  if (!_isValidForPotential()) return;
  _nVar = _ctxt.getNVar();
  _ctxt.setNVar(_nVar);

  // Initialize the working quantities
  _covpp = 0.;
  _covGp.resize(3);
  _covGG.resize(9);
  _dF.resize(3);
  _uF.resize(3);
  _hF.resize(3);
  _Tr.resize(9);
  _trttr.resize(9);
}

CovPotential::CovPotential(const CovPotential& r)
  : ACov(r)
  , _nVar(r._nVar)
  , _covRef(r._covRef)
  , _launchCalculations(r._launchCalculations)
  , _flagGradient(r._flagGradient)
  , _covpp(r._covpp)
  , _covGp(r._covGp)
  , _covGG(r._covGG)
  , _dF(r._dF)
  , _uF(r._uF)
  , _hF(r._hF)
  , _Tr(r._Tr)
  , _trttr(r._trttr)
{
}

CovPotential::~CovPotential()
{
}

bool CovPotential::_isValidForPotential() const
{
  auto ndim = _covRef.getNDim();
  auto nvar = _covRef.getNVar();
  if (ndim > 3)
  {
    messerr("This class is limited to 1-D or 2-D");
    return false;
  }
  if (nvar != 1)
  {
    messerr("This class is limited to Monovariate case");
    return false;
  }
  return true;
}

void CovPotential::_optimizationSetTarget(SpacePoint& pt) const
{
  DECLARE_UNUSED(pt)
}

/**
 * @brief According to the variable rank, call covariance between the variable and its derivatives
 *
 * @param p1 First point for covariance calculation
 * @param p2 Second point for covariance calculation
 * @param ivar Rank of the first variable (see remarks)
 * @param jvar Rank for the second variable (see remarks)
 * @param mode CovCalcMode structure
 * @return double
 *
 * @remark The use of this function is limited to the Monovariate case: therefore
 * @remark the arguments 'ivar' and 'jvar' cannot reference the variable rank.
 * @remark Instead the following convention applies for 'ivar' and 'jvar':
 * @remark - when =0, it refers to the variable itself
 * @remark - when >10, it applies to the gradient component (10 alongX, 11 along Y, ...)
 */
double CovPotential::_eval(const SpacePoint& p1,
                           const SpacePoint& p2,
                           Id ivar,
                           Id jvar,
                           const CovCalcMode* mode) const
{
  if (_launchCalculations)
    _evalZAndGradients(p1, p2);
  _launchCalculations = false;

  if (ivar == 0)
  {
    if (jvar == 0)
      return _covRef.evalCov(p1, p2, ivar, jvar, mode);

    int jdim = jvar - 1;
    return _covGp[jdim];
  }
  int idim = ivar - 1;
  if (jvar == 0)
    return _covGp[idim];
  int jdim = jvar - 1;
  return _covGG[3 * idim + jdim];
}

String CovPotential::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  sstr << "Covariance for Potential Model" << std::endl;
  return sstr.str();
}

/**
 * Calculate the square of the transformation matrix which transforms
 * a vector into its isotropic equivalent
 */
void CovPotential::_calculateTrTtr() const
{
  auto ndim = getNDim();
  _uF.fill(0.);
  _hF.fill(0.);
  _Tr.fill(0.);
  _trttr.fill(0.);

  // Matrix Tr = diag(coeffs) . R
  const MatrixSquare& mat    = _covRef.getAniso().getRotation().getMatrixDirect();
  const VectorDouble& scales = _covRef.getCorAniso()->getScales();

  for (Id i = 0; i < 3; i++)
    for (Id j = 0; j < 3; j++)
    {
      if (i >= ndim || j >= ndim) continue;
      TR(i, j) = mat.getValue(i, j) / scales[i];
    }

  // Calculate the t(Tr) %*% Tr matrix

  Id ecr = 0;
  for (Id i = 0; i < 3; i++)
    for (Id j = 0; j < 3; j++)
    {
      double prod = 0.;
      for (Id k = 0; k < 3; k++) prod += TR(k, i) * TR(k, j);
      _trttr[ecr++] = prod;
    }

  // Rotate the Vector

  _hF[0] = _Tr[0] * _dF[0] + _Tr[1] * _dF[1] + _Tr[2] * _dF[2];
  _hF[1] = _Tr[3] * _dF[0] + _Tr[4] * _dF[1] + _Tr[5] * _dF[2];
  _hF[2] = _Tr[6] * _dF[0] + _Tr[7] * _dF[1] + _Tr[8] * _dF[2];

  _uF[0] = _Tr[0] * _hF[0] + _Tr[3] * _hF[1] + _Tr[6] * _hF[2];
  _uF[1] = _Tr[1] * _hF[0] + _Tr[4] * _hF[1] + _Tr[7] * _hF[2];
  _uF[2] = _Tr[2] * _hF[0] + _Tr[5] * _hF[1] + _Tr[8] * _hF[2];
}

/**
 * Evaluates the covariance and gradient components
 * This function is restricted to the monovariate case
 * This function is limited to the only functions for which the
 * covariance(Point-Gradient) and covariance(Gradient-Gradient)
 * have been coded: i.e. Cubic or Gaussian, in addition to the nugget effect.
 *
 * If 'flag_grad' == 0, then the output array 'covGG' is not filled.
 *
 * @param p1     First point of the Increment
 * @param p2     Second point of the increment
 */

void CovPotential::_evalZAndGradients(const SpacePoint& p1,
                                      const SpacePoint& p2) const
{
  _covGp.fill(0);
  _covGG.fill(0);

  // Calculate the isotropic distance
  double hh = getSpace()->getDistance(p1, p2, _covRef.getAniso());

  //  Calculate the covariance
  double covar = _covRef.getSill(0, 0) * _covRef.getCorFunc()->evalCorFunc(hh);
  _covpp += covar;
  if (_covRef.getCorFunc()->getType() == ECov::NUGGET) return;

  // Calculate distance and plunge into a 3-D vector
  VectorDouble d1 = VH::subtract(p1.getCoords(), p2.getCoords());
  for (Id i = 0; i < 3; i++)
    _dF[i] = (i < static_cast<Id>(d1.size())) ? d1[i] : 0.;

  _calculateTrTtr();
  double dcovsr = _covRef.getSill(0, 0) * _covRef.getCorFunc()->evalCovDerivative(1, hh);

  //  Case where distance is null
  if (hh < EPSGRAD)
  {
    if (_flagGradient)
    {
      for (Id i = 0; i < 9; i++)
        _covGG[i] -= dcovsr * _trttr[i];
    }
  }
  else
  {
    //  Calculate covariance between point and gradient
    for (Id i = 0; i < 3; i++)
    {
      _covGp[i] += _uF[i] * dcovsr;
    }

    //  Calculate the covariance between gradient and gradient
    if (_flagGradient)
    {
      double d2cov = _covRef.getSill(0, 0) * _covRef.getCorFunc()->evalCovDerivative(2, hh);
      double a     = (dcovsr - d2cov) / (hh * hh);

      if (_covRef.getAniso().isIsotropic())
      {

        //  Isotropic case

        double b = dcovsr * _trttr[0];
        for (Id i = 0, ecr = 0; i < 3; i++)
          for (Id j = 0; j < 3; j++, ecr++)
          {
            _covGG[ecr] += a * _uF[i] * _uF[j];
            if (i == j) _covGG[ecr] -= b;
          }
      }
      else
      {

        //  Anisotropic case

        for (Id i = 0, ecr = 0; i < 3; i++)
          for (Id j = 0; j < 3; j++, ecr++)
          {
            _covGG[ecr] += a * _uF[i] * _uF[j] - dcovsr * _trttr[ecr];
          }
      }
    }
  }
}

} // namespace gstlrn