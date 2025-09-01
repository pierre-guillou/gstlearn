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
#include "Covariances/CovGradientGeneric.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

#include <cmath>

namespace gstlrn
{
CovGradientGeneric::CovGradientGeneric(const ACov& cova, double ballradius)
  : ACov()
  , _ballRadius(ballradius)
  , _covRef(cova)
{
  setContext(cova.getContext());
  if (!_isValidForGradient()) return;

  // Consider the variable and its gradient(s)
  Id nVar = _ctxt.getNVar() + static_cast<Id>(_ctxt.getNDim());
  _ctxt.setNVar(nVar);
}

CovGradientGeneric::CovGradientGeneric(const CovGradientGeneric& r)
  : ACov(r)
  , _ballRadius(r._ballRadius)
  , _covRef(r._covRef)
{
}

CovGradientGeneric::~CovGradientGeneric()
{
}

bool CovGradientGeneric::_isValidForGradient() const
{
  auto ndim = _covRef.getNDim();
  auto nvar = _covRef.getNVar();
  if (ndim > 2)
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

void CovGradientGeneric::_optimizationSetTarget(SpacePoint& pt) const
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
 * @remark This use of this function is limited to the Monovariate case.
 * @remark The argument 'ivar' (resp. 'jvar') gives the variable rank (if equal to 0)
 * Otherwise it gives the space index derivative (=idim-1)
 */
double CovGradientGeneric::_eval(const SpacePoint& p1,
                                 const SpacePoint& p2,
                                 Id ivar,
                                 Id jvar,
                                 const CovCalcMode* mode) const
{
  if (ivar == 0 && jvar == 0)
    return _covRef.evalCov(p1, p2, ivar, jvar, mode);

  Id idim = ivar - 1;
  Id jdim = jvar - 1;
  if (ivar == 0)
    return -_evalZGradientNumeric(p1, p2, jdim, _ballRadius, mode);
  if (jvar == 0)
    return +_evalZGradientNumeric(p1, p2, idim, _ballRadius, mode);
  if (jdim == idim)
    return _evalGradientGradientNumeric(p1, p2, idim, jdim, _ballRadius, mode);
  return -_evalGradientGradientNumeric(p1, p2, idim, jdim, _ballRadius, mode);
}

String CovGradientGeneric::toString(const AStringFormat* strfmt) const
{
  DECLARE_UNUSED(strfmt);
  std::stringstream sstr;

  sstr << "Covariance for Variable and its Gradients" << std::endl;
  sstr << "Derivation distance" << _ballRadius << std::endl;
  return sstr.str();
}

double CovGradientGeneric::_evalZGradientNumeric(const SpacePoint& p1,
                                                 const SpacePoint& p2,
                                                 Id idim,
                                                 double radius,
                                                 const CovCalcMode* mode) const
{
  SpacePoint paux;
  auto ndim = getContext().getNDim();
  VectorDouble vec(ndim, 0);

  vec[idim] = +radius / 2.;
  paux      = p2;
  paux.move(vec);
  double covp0 = _covRef.evalCov(p1, paux, 0, 0, mode);
  vec[idim]    = -radius / 2.;
  paux         = p2;
  paux.move(vec);
  double covm0 = _covRef.evalCov(p1, paux, 0, 0, mode);

  double cov = (covm0 - covp0) / radius;
  return (cov);
}

double CovGradientGeneric::_evalGradientGradientNumeric(const SpacePoint& p1,
                                                        const SpacePoint& p2,
                                                        Id idim,
                                                        Id jdim,
                                                        double radius,
                                                        const CovCalcMode* mode) const
{
  SpacePoint paux;
  auto ndim = getContext().getNDim();
  VectorDouble vec(ndim, 0);

  double cov;
  if (idim != jdim)
  {
    vec[idim] = -radius / 2.;
    vec[jdim] = +radius / 2.;
    paux      = p2;
    paux.move(vec);
    double covmp = _covRef.evalCov(p1, paux, 0, 0, mode);
    vec[idim]    = -radius / 2.;
    vec[jdim]    = -radius / 2.;
    paux         = p2;
    paux.move(vec);
    double covmm = _covRef.evalCov(p1, paux, 0, 0, mode);
    vec[idim]    = +radius / 2.;
    vec[jdim]    = -radius / 2.;
    paux         = p2;
    paux.move(vec);
    double covpm = _covRef.evalCov(p1, paux, 0, 0, mode);
    vec[idim]    = +radius / 2.;
    vec[jdim]    = +radius / 2.;
    paux         = p2;
    paux.move(vec);
    double covpp = _covRef.evalCov(p1, paux, 0, 0, mode);

    cov = (covmm + covpp - covmp - covpm) / (radius * radius);
  }
  else
  {
    vec[idim] = +radius;
    paux      = p2;
    paux.move(vec);
    double cov2m = _covRef.evalCov(p1, paux, 0, 0, mode);
    vec[idim]    = -radius;
    paux         = p2;
    paux.move(vec);
    double cov2p = _covRef.evalCov(p1, paux, 0, 0, mode);
    vec[idim]    = 0;
    paux         = p2;
    paux.move(vec);
    double cov00 = _covRef.evalCov(p1, paux, 0, 0, mode);

    cov = -2. * (cov2p - 2. * cov00 + cov2m) / (radius * radius);
  }
  return (cov);
}

} // namespace gstlrn