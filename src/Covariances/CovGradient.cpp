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
#include "Covariances/CovGradient.hpp"
#include "Covariances/CovAniso.hpp"
#include "Covariances/CovCalcMode.hpp"
#include "Covariances/CovContext.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpaceComposite.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"

#include <cmath>
#include <vector>

namespace gstlrn
{
CovGradient::CovGradient(const CovAniso& cova)
  : ACov()
  , _nVar(0) // to be fixed
  , _covRef(cova)
{
  setContext(_covRef.getContext());
  _ctxt.setNVar(_nVar);
}

CovGradient::CovGradient(const CovGradient& r)
  : ACov(r)
  , _nVar(r._nVar)
  , _covRef(r._covRef)
{
}

CovGradient& CovGradient::operator=(const CovGradient& r)
{
  if (this != &r)
  {
    ACov::operator=(r);
    _nVar = r._nVar;
    // The member _covRef cannot be updated by this method
  }
  return *this;
}

CovGradient::~CovGradient()
{
}

void CovGradient::_optimizationSetTarget(SpacePoint& pt) const
{
  DECLARE_UNUSED(pt)
}

void CovGradient::_optimizationPreProcess(int mode, const std::vector<SpacePoint>& ps) const
{
  DECLARE_UNUSED(mode)
  DECLARE_UNUSED(ps)
}

void CovGradient::_optimizationPostProcess() const
{
}

double CovGradient::_eval(const SpacePoint& p1,
                          const SpacePoint& p2,
                          int ivar,
                          int jvar,
                          const CovCalcMode* mode) const
{
  return _covRef._eval(p1, p2, ivar, jvar, mode);
}

} // namespace gstlrn