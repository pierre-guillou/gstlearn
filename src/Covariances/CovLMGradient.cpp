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
#include "Covariances/CovLMGradient.hpp"
#include "Covariances/CovAnisoList.hpp"
#include "Covariances/CovGradientFunctional.hpp"
#include "Space/ASpace.hpp"
#include "Space/SpacePoint.hpp"

namespace gstlrn
{
CovLMGradient::CovLMGradient(const CovContext& ctxt)
  : CovAnisoList(ctxt)
{
  setOptimEnabled(false);
}

CovLMGradient::CovLMGradient(const CovLMGradient& r)
  : CovAnisoList(r)
{
  setOptimEnabled(false);
}

CovLMGradient::CovLMGradient(const CovAnisoList& r)
  : CovAnisoList(r.getContext())
{
  setOptimEnabled(false);
  for (auto icov = r.getNCov() - 1; icov >= 0; icov--)
  {
    const CovAniso* cov = r.getCovAniso(icov);
    if (!cov->hasCovDerivative())
    {
      messerr("The covariance %s is not compatible with Gradients",
              cov->getCovName().c_str());
    }
    else
    {
      CovGradientFunctional newcov(*cov);
      addCov(newcov);
    }
  }
  for (auto& e: _covs)
  {
    static_cast<CovAniso*>(e.get())->setOptimEnabled(false);
  }
}

CovLMGradient& CovLMGradient::operator=(const CovLMGradient& r)
{
  if (this != &r)
  {
    CovAnisoList::operator=(r);
  }
  setOptimEnabled(false);
  return *this;
}

CovLMGradient::~CovLMGradient()
{
  /// TODO : Delete pointers ?
}

void CovLMGradient::evalZAndGradients(const SpacePoint& p1,
                                      const SpacePoint& p2,
                                      double& covVal,
                                      VectorDouble& covGp,
                                      VectorDouble& covGG,
                                      const CovCalcMode* mode,
                                      bool flagGrad) const
{
  covVal = 0.;
  covGp.fill(0.);
  if (flagGrad)
    covGG.fill(0.);

  for (size_t i = 0, n = getNCov(); i < n; i++)
  {
    auto* covloc = dynamic_cast<ACovGradient*>(_covs[i].get());
    if (covloc != nullptr)
      covloc->evalZAndGradients(p1, p2, covVal, covGp, covGG, mode, flagGrad);
  }
}

void CovLMGradient::evalZAndGradients(const VectorDouble& vec,
                                      double& covVal,
                                      VectorDouble& covGp,
                                      VectorDouble& covGG,
                                      const CovCalcMode* mode,
                                      bool flagGrad) const
{
  /// TODO : Not true whatever the space
  SpacePoint p1(getSpace()->getOrigin(), -1);
  SpacePoint p2(getSpace()->getOrigin(), -1);
  p2.move(vec);

  evalZAndGradients(p1, p2, covVal, covGp, covGG, mode, flagGrad);
}

void CovLMGradient::addCov(const CovBase& cov)
{
  // TODO This should be checked for some cases of Gradient (probably non numerical)
  const auto* covgrad = dynamic_cast<const ACovGradient*>(&cov);
  if (covgrad == nullptr)
  {
    messerr("This covariance cannot be added");
    return;
  }
  cov.setOptimEnabled(false);
  CovAnisoList::addCov(cov);
}

} // namespace gstlrn