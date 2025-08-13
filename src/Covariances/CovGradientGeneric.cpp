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
  , _nVar(0)
  , _ballRadius(ballradius)
  , _covRef(cova)
{
  setContext(cova.getContext());
  if (!_isValidForGradient()) return;
  _nVar = _ctxt.getNVar() + _ctxt.getNDim(); // Consider the variable and its gradient(s)
  _ctxt.setNVar(_nVar);
}

CovGradientGeneric::CovGradientGeneric(const CovGradientGeneric& r)
  : ACov(r)
  , _nVar(r._nVar)
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

// void CovGradientGeneric::_optimizationPreProcess(Id mode, const std::vector<SpacePoint>& ps) const
// {
//   DECLARE_UNUSED(mode)
//   DECLARE_UNUSED(ps)
// }

// void CovGradientGeneric::_optimizationPostProcess() const
// {
// }

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
    return -_covRef.evalZGNumeric(p1, p2, 0, 0, jdim, _ballRadius, mode);
  if (jvar == 0)
    return +_covRef.evalZGNumeric(p1, p2, 0, 0, idim, _ballRadius, mode);
  if (jdim == idim)
    return _covRef.evalGGNumeric(p1, p2, 0, 0, idim, jdim, _ballRadius, mode);
  return -_covRef.evalGGNumeric(p1, p2, 0, 0, idim, jdim, _ballRadius, mode);
}

} // namespace gstlrn