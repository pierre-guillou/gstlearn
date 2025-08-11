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

#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CorAniso.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Space/SpacePoint.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include <vector>

namespace gstlrn
{
class ACov;
class CovAniso;
/**
 * \brief
 * This class describes the Gneiting correlation function.
 *
 */
class GSTLEARN_EXPORT CovGradient: public ACov
{
public:
  CovGradient(const CovAniso& cova);
  CovGradient(const CovGradient& r);
  CovGradient& operator=(const CovGradient& r);
  virtual ~CovGradient();
  IMPLEMENT_CLONING(CovGradient)

  bool isConsistent(const ASpace* space) const override
  {
    DECLARE_UNUSED(space)
    return true;
  }

  /// ACov Interface
  int getNVar() const override { return _nVar; }

protected:
  double _eval(const SpacePoint& p1,
               const SpacePoint& p2,
               int ivar                = 0,
               int jvar                = 0,
               const CovCalcMode* mode = nullptr) const override;
  void _optimizationSetTarget(SpacePoint& pt) const override;

private:
  void _optimizationPreProcess(int mode, const std::vector<SpacePoint>& ps) const override;
  void _optimizationPostProcess() const override;

private:
  int _nVar; // TODO should be number of variables and gradients
  const CovAniso& _covRef;
};
} // namespace gstlrn
