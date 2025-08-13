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

#include "Covariances/ACov.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class ACov;
class CovAniso;
/**
 * \brief
 * This class describes the Covariance to be used when processing Data
 * and its gradient components.
 * This covariance is based on the initial covariance of the Data and derives
 * the simple and cross covariances of its gradient components.
 * It uses Numerical Derivation and therefore is suitable whatever the type
 * of covariance used for the Data variable.
 */
class GSTLEARN_EXPORT CovGradientGeneric: public ACov
{
public:
  CovGradientGeneric(const ACov& cova, double ballradius = 1);
  CovGradientGeneric(const CovGradientGeneric& r);
  CovGradientGeneric& operator=(const CovGradientGeneric& r) = delete;
  virtual ~CovGradientGeneric();

  /// ICloneable Interface
  IMPLEMENT_CLONING(CovGradientGeneric)

  bool isConsistent(const ASpace* space) const override
  {
    DECLARE_UNUSED(space)
    return true;
  }

  /// ACov Interface
  Id getNVar() const override { return _nVar; }

protected:
  double _eval(const SpacePoint& p1,
               const SpacePoint& p2,
               Id ivar                 = 0,
               Id jvar                 = 0,
               const CovCalcMode* mode = nullptr) const override;
  void _optimizationSetTarget(SpacePoint& pt) const override;

private:
  // void _optimizationPreProcess(Id mode, const std::vector<SpacePoint>& ps) const override;
  // void _optimizationPostProcess() const override;
  bool _isValidForGradient() const;

private:
  Id _nVar;
  double _ballRadius;
  const ACov& _covRef;
};

} // namespace gstlrn
