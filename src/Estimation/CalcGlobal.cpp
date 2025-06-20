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
#include "Covariances/CovCalcMode.hpp"
#include "Enum/ECalcMember.hpp"
#include "geoslib_old_f.h"

#include "Estimation/CalcGlobal.hpp"
#include "Db/DbGrid.hpp"
#include "Db/Db.hpp"
#include "Estimation/KrigingSystem.hpp"
#include "Neigh/NeighUnique.hpp"
#include "Model/Model.hpp"

#include <math.h>

CalcGlobal::CalcGlobal(int ivar0, bool verbose)
  : ACalcInterpolator()
  , _flagArithmetic(false)
  , _flagKriging(false)
  , _ivar0(ivar0)
  , _verbose(verbose)
  , _modelLocal(nullptr)
{
}

CalcGlobal::~CalcGlobal()
{
}

bool CalcGlobal::_check()
{
  if (! ACalcInterpolator::_check()) return false;

  if (! hasDbin()) return false;
  if (! hasDbout()) return false;
  if (! hasModel()) return false;

  _modelLocal = dynamic_cast<Model*>(getModel());
  if (_modelLocal == nullptr)
  {
    messerr("This method requires the model to be a 'Model' (not a ModelGeneric)");
    return false;
  }

  if (_flagArithmetic)
  {
    if (! getDbout()->isGrid())
    {
      messerr("'dbout'  must be a grid for Arithmetic Global estimation");
      return false;
    }
  }

  if (_ivar0 < 0 || _ivar0 >= getDbin()->getNLoc(ELoc::Z))
  {
    messerr("The target variable (%d) must lie between 1 and the number of variables (%d)",
            _ivar0 + 1, getDbin()->getNLoc(ELoc::Z));
    return false;
  }

  return true;
}

bool CalcGlobal::_preprocess()
{
  return ACalcInterpolator::_preprocess();
}

bool CalcGlobal::_postprocess()
{
  _cleanVariableDb(2);

  return true;
}

void CalcGlobal::_rollback()
{
  _cleanVariableDb(1);
}

/****************************************************************************/
/*!
 **  Standard Kriging
 **
 ** \return  Error return code
 **
 *****************************************************************************/
bool CalcGlobal::_run()
{
  if (_flagArithmetic)
  {
    if (_globalArithmetic()) return false;
  }

  if (_flagKriging)
  {
    if (_globalKriging()) return false;
  }
  return true;
}

int CalcGlobal::_globalKriging()
{
  VectorDouble rhsCum;
  Db* dbin            = getDbin();
  Db* dbout           = getDbout();
  int nvar            = _modelLocal->getNVar();
  int ng              = 0;
  VectorDouble wgt;

  KrigOpt krigopt;
  MatrixSymmetric Sigma;
  MatrixDense X;

  // Get the Covariance between data (Unique Neighborhood)
  CovCalcMode mode            = CovCalcMode(ECalcMember::LHS);
  VectorVectorInt sampleRanks = dbin->getSampleRanks({_ivar0});
  VectorDouble Z              = dbin->getValuesByRanks(sampleRanks,
                                                       _modelLocal->getMeans(), 
                                                       !_modelLocal->hasDrift());
  if (_modelLocal->evalCovMatSymInPlaceFromIdx(Sigma, dbin, sampleRanks, &mode, false)) return 1;
  if (_modelLocal->evalDriftMatByRanksInPlace(X, dbin, sampleRanks, ECalcMember::LHS)) return 1;

  KrigingAlgebra algebra;
  algebra.resetNewData();
  algebra.setData(&Z, &sampleRanks, &_modelLocal->getMeans());
  algebra.setLHS(&Sigma, &X);

  // Prepare the cumulative matrices
  MatrixDense Sigma0Cum(Sigma.getNRows(), 1);
  MatrixDense X0Cum(1, X.getNCols());
  MatrixDense Sigma0;
  MatrixDense X0;
  MatrixSymmetric Sigma00;

  /* Loop on the targets to be processed */
  for (int iech = 0, nech = dbout->getNSample(); iech < nech; iech++)
  {
    mes_process("Kriging sample", dbout->getNSample(), iech);
    if (!dbout->isActive(iech)) continue;

    if (_modelLocal->evalCovMatRHSInPlaceFromIdx(Sigma0, dbin, dbout, sampleRanks, iech, krigopt, false)) return 1;
    if (_modelLocal->evalDriftMatByTargetInPlace(X0, dbout, iech, krigopt)) return 1;

    // Cumulate the R.H.S.
    Sigma0Cum.addMatInPlace(Sigma0);
    X0Cum.addMatInPlace(X0);
    ng++;
  }

  // Normalize the cumulative R.H.S.
  double oneOverNG = 1. / (double)ng;
  Sigma0Cum.prodScalar(oneOverNG);
  X0Cum.prodScalar(oneOverNG);
  algebra.setRHS(&Sigma0Cum, &X0Cum);

  if (_modelLocal->evalCovMat0InPlace(Sigma00, dbout, 0)) return 1;
  algebra.setVariance(&Sigma00);

  double estim = algebra.getEstimation()[0];
  double stdv  = algebra.getStdv()[0];
  // The previous term corresponds to the standard deviation calculated
  // with a punctual target. Therefore the corresponding variance
  // must be corrected (C00 -> Cvv) to pass to a correct Territory variance 
  // of estimation
  double c00   = Sigma00.getValue(0,0);

  /* Preliminary checks */

  int ntot = dbin->getNSample(false);
  int np   = dbin->getNSample(true);
  double cell = 1.;
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(dbout);
  if (dbgrid != nullptr) cell = dbgrid->getCellSize();
  double surface = ng * cell;

  /* Average covariance over the territory */

  double cvv = _modelLocal->evalAverageDbToDb(dbout, dbout, _ivar0, _ivar0,
                                        dbin->getExtensionDiagonal() / 1.e3, 0);

  /* Perform the estimation */

  double cvvgeo = stdv * stdv - c00 + cvv;

  double stdgeo = (cvvgeo > 0) ? sqrt(cvvgeo) : 0.;
  double cvgeo = (isZero(estim) || FFFF(estim)) ? TEST : stdgeo / estim;

  /* Store the results in the output Global_Result struture */

  _gRes.ntot = ntot;
  _gRes.np = np;
  _gRes.ng = ng;
  _gRes.surface = surface;
  _gRes.zest = estim;
  _gRes.sse = stdgeo;
  _gRes.cvgeo = cvgeo;
  _gRes.cvv = cvv;
  _gRes.weights = wgt;

  /* Printout */

  if (_verbose)
  {
    mestitle(1,"Global estimation kriging");
    message("Total number of data             = %d\n", ntot);
    message("Number of active data            = %d\n", np);
    message("Number of variables              = %d\n", nvar);
    message("Cvv                              = %lf\n", cvv);
    if (FFFF(estim))
      message("Estimation by kriging            = NA\n");
    else
      message("Estimation by kriging            = %lf\n", estim);
    message("Estimation St. Dev. of the mean  = %lf\n", stdgeo);
    if (FFFF(cvgeo))
      message("CVgeo                            = NA\n");
    else
      message("CVgeo                            = %lf\n", cvgeo);
    message("Surface                          = %lf\n", surface);
    if (FFFF(estim))
      message("Q (Estimation * Surface)         = NA\n");
    else
      message("Q (Estimation * Surface)         = %lf\n", estim * surface);
    message("\n");
  }
  _modelLocal->optimizationPostProcess();
  return 0;
}

int CalcGlobal::_globalArithmetic()
{
  DbGrid* dbgrid = dynamic_cast<DbGrid*>(getDbout());
  int ntot = getDbin()->getNSample(false);
  int np = getDbin()->getNSample(true);
  int ng = dbgrid->getNSample(true);
  double surface = ng * dbgrid->getCellSize();

  /* Average covariance over the data */

  double cxx =
    _modelLocal->evalAverageDbToDb(getDbin(), getDbin(), _ivar0, _ivar0, 0., 0);

  /* Average covariance between the data and the territory */

  double cxv =
    _modelLocal->evalAverageDbToDb(getDbin(), dbgrid, _ivar0, _ivar0, 0., 0);

  /* Average covariance over the territory */

  double cvv = _modelLocal->evalAverageDbToDb(dbgrid, dbgrid, _ivar0, _ivar0,
                                        dbgrid->getExtensionDiagonal() / 1.e3, 0);

  /* Calculating basic statistics */

  int iatt = getDbin()->getUIDByLocator(ELoc::Z, _ivar0);
  double wtot;
  double ave;
  double var;
  double mini;
  double maxi;
  db_monostat(getDbin(), iatt, &wtot, &ave, &var, &mini, &maxi);

  /* Filling the resulting structure */

  double sse = cvv - 2. * cxv + cxx;
  sse = (sse > 0) ? sqrt(sse) : 0.;
  double cvsam = (! isZero(ave)) ? sqrt(var) / ave : TEST;
  double cviid = cvsam / sqrt(np);
  double cvgeo = (! isZero(ave)) ? sse / ave : TEST;

  /* Filling the output structure */

  _gRes.ntot = ntot;
  _gRes.np = np;
  _gRes.ng = ng;
  _gRes.surface = surface;
  _gRes.zest = ave;
  _gRes.sse  = sse;
  _gRes.cvgeo = cvgeo;
  _gRes.cvv = cvv;
  _gRes.weights.resize(np, 1./np);

  if (_verbose)
  {
    mestitle(1,"Global estimation by arithmetic average");
    message("Total number of data             = %d\n", ntot);
    message("Number of active data            = %d\n", np);
    message("Sample variance                  = %lf\n", var);
    message("CVsample                         = %lf\n", cvsam);
    message("CViid                            = %lf\n", cviid);
    message("Cxx                              = %lf\n", cxx);
    message("Cxv                              = %lf\n", cxv);
    message("Cvv                              = %lf\n", cvv);
    if (FFFF(ave))
      message("Estimation by arithmetic average = NA\n");
    else
      message("Estimation by arithmetic average = %lf\n", ave);
    message("Estimation St. dev. of the mean  = %lf\n", sse);
    if (FFFF(cvgeo))
      message("CVgeo                            = NA\n");
    else
      message("CVgeo                            = %lf\n", cvgeo);
    message("Surface                          = %lf\n", surface);
    if (FFFF(ave))
      message("Q (Estimation * Surface)         = NA\n");
    else
      message("Q (Estimation * Surface)         = %lf\n", ave * surface);
    message("\n");
  }

  return 0;
}

Global_Result global_arithmetic(Db *dbin,
                                DbGrid *dbgrid,
                                ModelGeneric *model,
                                int ivar0,
                                bool verbose)
{
  Global_Result gres;
  CalcGlobal global(ivar0, verbose);
  global.setDbin(dbin);
  global.setDbout(dbgrid);
  global.setModel(model);
  global.setFlagArithmetic(true);

  if (global.run())
    gres = global.getGRes();
  return gres;
}

Global_Result global_kriging(Db *dbin,
                            Db *dbout,
                            ModelGeneric *model,
                            int ivar0,
                            bool verbose)
{
  Global_Result gres;
  CalcGlobal global(ivar0, verbose);
  global.setDbin(dbin);
  global.setDbout(dbout);
  global.setModel(model);
  global.setFlagKriging(true);

  if (global.run())
    gres = global.getGRes();
  return gres;
}

