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
#include "Gibbs/GibbsUPropMono.hpp"
#include "Model/Model.hpp"
#include "Db/Db.hpp"
#include "Basic/Law.hpp"
#include "Basic/OptDbg.hpp"
#include "Model/CovInternal.hpp"

#include <math.h>

GibbsUPropMono::GibbsUPropMono()
  : GibbsMultiMono()
  , _rval(0.5)
  , _eps(EPSILON3)
{
}

GibbsUPropMono::GibbsUPropMono(Db* db, const std::vector<Model *>& models, double rho)
  : GibbsMultiMono(db, models, rho)
  , _rval(0.5)
  , _eps(EPSILON3)
{
}

GibbsUPropMono::GibbsUPropMono(const GibbsUPropMono &r)
  : GibbsMultiMono(r)
  , _rval(0.5)
  , _eps(r._eps)
{
}

GibbsUPropMono& GibbsUPropMono::operator=(const GibbsUPropMono &r)
{
  if (this != &r)
  {
    AGibbs::operator=(r);
    _rval = r._rval;
    _eps = r._eps;
  }
  return *this;
}

GibbsUPropMono::~GibbsUPropMono()
{
}

/****************************************************************************/
/*!
**  Establish the covariance matrix for Gibbs
**
** \return  Error returned code
**
** \param[in]  verbose      Verbose flag
** \param[in]  verboseTimer True to show elapse times
**
*****************************************************************************/
int GibbsUPropMono::covmatAlloc(bool verbose, bool /*verboseTimer*/)
{
  if (verbose) mestitle(1,"Gibbs using Unique Neighborhood in Propagative case");

  // Initialize the statistics (optional)

  _statsInit();

  return 0;
}

/****************************************************************************/
/*!
**  Perform one update of the Gibbs sampler (Propagative algorithm)
**
** \param[in]  y           Gaussian vector
** \param[in]  isimu       Rank of the simulation
** \param[in]  ipgs        Rank of the GS (should be 0)
** \param[in]  iter        Rank of the iteration
**
*****************************************************************************/
void GibbsUPropMono::update(VectorVectorDouble& y,
                            int isimu,
                            int ipgs,
                            int iter)
{
  CovCalcMode mode;

  /* Initializations */

  Db* db = getDb();
  Model* model = getModels(0);
  int nact  = _getSampleRankNumber();
  int ndim  = model->getNDim();
  int icase = getRank(ipgs,0);

  double eps  = getEps();
  double r    = getRval();
  double sqr  = sqrt(1. - r * r);

  /* Core allocation */

  VectorDouble d1(ndim);
  VectorBool img(nact * nact);

  /* Print the title */

  if (OptDbg::query(EDbg::CONVERGE))
    mestitle(1,"Iterative Conditional Expectation (Simu:%d)",isimu+1);

  /* Loop on the samples */

  for (int iact = 0; iact < nact; iact++)
  {
    int iech = getSampleRank(iact);

    /* Covariance vector between the current datum and the other samples */

    double sigval;
    for (int idim = 0; idim < ndim; idim++)
      d1[idim] = 0.;
    if (model->getCovAnisoList()->isNoStat())
    {
      CovInternal covint(1, iech, 1, iech, ndim, db, db);
      sigval = model->evaluateOneGeneric(&covint, d1);
    }
    else
    {
      sigval = model->evaluateOneGeneric(nullptr, d1);
    }
    if (sigval <= 0) continue;
    sigval = sqrt(sigval);
    double delta = (r - 1.) * y[icase][iact] + sigval * sqr * law_gaussian();

    /* Update the gaussian vector */

    for (int jact = 0; jact < nact; jact++)
    {
      if (iter > 0 && ! img[nact * iact + jact]) continue;
      int jech = getSampleRank(jact);

      double sigloc;
      for (int idim = 0; idim < ndim; idim++)
        d1[idim] = db->getCoordinate(iech, idim) - db->getCoordinate(jech, idim);
      if (model->getCovAnisoList()->isNoStat())
      {
        CovInternal covint(1, iech, 1, jech, ndim, db, db);
        sigloc = model->evaluateOneGeneric(&covint, d1);
      }
      else
      {
        sigloc = model->evaluateOneGeneric(nullptr, d1);
      }

      bool flag_affect = (ABS(sigloc) > sigval * eps);
      if (iter <= 0) img[nact * iact + jact] = flag_affect;
      if (flag_affect) y[icase][jact] += delta * sigloc / sigval;
    }
  }

  // Update statistics (optional)

  _updateStats(y, ipgs, iter);
}
