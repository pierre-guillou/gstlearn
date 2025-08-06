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

#include "Basic/ICloneable.hpp"
#include "Covariances/ACovGradient.hpp"
#include "Covariances/CovContext.hpp"

namespace gstlrn
{
class Rotation;

/**
 * Class dedicated to manipulating a variables and its derivatives.
 * This feature is limited to the monovariate case
 */
class GSTLEARN_EXPORT CovGradientNumerical: public ACovGradient
{
public:
  CovGradientNumerical(const ECov& type, double ballRadius, const CovContext& ctxt);
  CovGradientNumerical(const CovGradientNumerical& r);
  CovGradientNumerical& operator=(const CovGradientNumerical& r);
  virtual ~CovGradientNumerical();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovGradientNumerical)

  virtual double eval0(Id ivar = 0,
                       Id jvar = 0,
                       const CovCalcMode* mode = nullptr) const override;


  double getBallRadius() const override { return _ballRadius; }

  void evalZAndGradients(const SpacePoint& p1,
                         const SpacePoint& p2,
                         double& covVal,
                         VectorDouble& covGp,
                         VectorDouble& covGG,
                         const CovCalcMode* mode = nullptr,
                         bool flagGrad = false) const override;

protected:

virtual double _eval(const SpacePoint& p1,
                     const SpacePoint& p2,
                     Id ivar = 0,
                     Id jvar = 0,
                     const CovCalcMode* mode = nullptr) const override;
private:
  double _evalZZ(Id ivar,
                 Id jvar,
                 const SpacePoint& p1,
                 const SpacePoint& p2,
                 const CovCalcMode* mode = nullptr) const;
  double _evalZGrad(Id ivar,
                    Id jvar,
                    Id idim,
                    const SpacePoint& p1,
                    const SpacePoint& p2,
                    const CovCalcMode* mode = nullptr) const;
  double _evalGradGrad(Id ivar,
                       Id jvar,
                       Id idim,
                       Id jdim,
                       const SpacePoint& p1,
                       const SpacePoint& p2,
                       const CovCalcMode* mode = nullptr) const;

private:
  double _ballRadius;   /*! Radius of the Ball for Numerical Gradient calculation */
};

}