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
 * This class describes the Covariance to be used when processing
 * Data, Gradients and Tangents in the Potential Method framework.
 * This covariance is based on the initial covariance of the Data and derives
 * all simple and cross covariances for using Gradients and/or Tangents.
 * It uses Functional Derivation and therefore is suitable for a limited
 * set of differentiable covariances.
 */
class GSTLEARN_EXPORT CovPotential: public ACov
{
public:
  CovPotential(const CovAniso& cova, bool flagGradient);
  CovPotential(const CovPotential& r);
  CovPotential& operator=(const CovPotential& r) = delete;
  virtual ~CovPotential();

  /// ICloneable Interface
  IMPLEMENT_CLONING(CovPotential)

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

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
  bool _isValidForPotential() const;
  void _checkPointHasChanged(const SpacePoint& p1,
                             const SpacePoint& p2) const;
  void _calculateTrTtr(const VectorDouble& d) const;
  void _evalZAndGradients(const SpacePoint& p1, const SpacePoint& p2) const;

private:
  Id _nVar;
  bool _flagGradient;
  const CovAniso& _covRef;

  mutable SpacePoint _p1Mem;
  mutable SpacePoint _p2Mem;
  // covpp  Covariance value
  mutable double _covpp;
  // covGp  Covariance <G[i](x0+x,y0+y,z0+z), P(x0,y0,z0)> (dim=3)
  mutable VectorDouble _covGp;
  // covGG  Covariance <G[i](x0+x,y0+y,z0+z), G[j](x0,y0,z0)> (dim=3)
  mutable VectorDouble _covGG;
  mutable VectorDouble _dF;
  mutable VectorDouble _uF;
  mutable VectorDouble _hF;
  mutable VectorDouble _Tr;
  mutable VectorDouble _trttr;
};

} // namespace gstlrn
