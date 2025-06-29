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

#include "Model/AModelFitSills.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Covariances/ACov.hpp"
#include "Covariances/CovContext.hpp"
#include "Drifts/DriftList.hpp"
#include "Basic/ListParams.hpp"
#include "Basic/ICloneable.hpp"
#include "Db/RankHandler.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Model/Constraints.hpp"
#include "Model/ModelOptimParam.hpp"

class Model;
class Db;

class DbGrid;
class CovCalcMode;
/**
 * \brief
 * Class containing the Model Information describing the formal Spatial (or Temporal) Characteristics
 * of the (set of) random variable(s) under study.
 *
 * The Model is essentially a container with two main contents:
 * - the **covariance** part: see ACov.hpp for more information
 * - the **drift** part: see DriftList.hpp for more information
 *
 * The additional member **CovContext** only serves in carrying the following information:
 * - the number of variables: if more than 1, the Model becomes multivariate
 * - the field extension: this information is needed to get a *stationary* version to any covariance
 * - the experimental mean vector and the variance-covariance matrix (used to calibrate the Model)
 */
class GSTLEARN_EXPORT ModelGeneric : public ICloneable
{
public:
  ModelGeneric(const CovContext& ctxt = CovContext());
  ModelGeneric(const ModelGeneric &r);
  ModelGeneric& operator= (const ModelGeneric &r);
  virtual ~ModelGeneric();

  //getters for member pointers
  const ACov*       getCov()             const { return  _cova;     }
  const CovContext* getContext()         const { return &_ctxt;     }
  const DriftList*  getDriftList()       const { return  _driftList;}

public:
  /// ICloneable interface
  IMPLEMENT_CLONING(ModelGeneric)

  ACov*       _getCovModify() { return _cova; }
  CovContext* _getContextModify() { return &_ctxt; }
  DriftList*  _getDriftListModify() { return _driftList; }
  
public:
  // Forwarding the methods from _cova
  FORWARD_METHOD(getCov, evalCovMat)
  FORWARD_METHOD(getCov, evalCovMatInPlace)
  FORWARD_METHOD(getCov, evalCovMatInPlaceFromIdx)
  FORWARD_METHOD(getCov, evalCovMatSym)
  FORWARD_METHOD(getCov, evalCovMatSymInPlace)
  FORWARD_METHOD(getCov, evalCovMatSymInPlaceFromIdx)
  FORWARD_METHOD(getCov, eval0Mat)
  FORWARD_METHOD(getCov, evalCovMat0)
  FORWARD_METHOD(getCov, evalCovMat0InPlace)
  FORWARD_METHOD(getCov, evalCovVecRHSInPlace)
  FORWARD_METHOD(getCov, evalCovMatOptimInPlace)
  FORWARD_METHOD(getCov, evalCovMatRHSInPlaceFromIdx)
  FORWARD_METHOD(getCov, evalCovMatSparse)
  FORWARD_METHOD(getCov, eval0)
  FORWARD_METHOD(getCov, evalCov)
  FORWARD_METHOD(getCov, evalNvarIpas)
  FORWARD_METHOD(getCov, evalNvarIpasIncr)
  FORWARD_METHOD(getCov, evalIvarNlag)
  FORWARD_METHOD(getCov, evalIvarIpas)
  FORWARD_METHOD(getCov, evalCvv)
  FORWARD_METHOD(getCov, evalCvvShift)
  FORWARD_METHOD(getCov, evalCvvM)
  FORWARD_METHOD(getCov, evalCxv)
  FORWARD_METHOD(getCov, evalCxvM)
  FORWARD_METHOD(getCov, evalPointToDb)
  FORWARD_METHOD(getCov, evalPointToDbAsSP)
  FORWARD_METHOD(getCov, evalAverageDbToDb,TEST)
  FORWARD_METHOD(getCov, evalAverageIncrToIncr,TEST)
  FORWARD_METHOD(getCov, evalAveragePointToDb,TEST)
  FORWARD_METHOD(getCov, samplingDensityVariance, TEST)
  FORWARD_METHOD(getCov, specificVolume, TEST)
  FORWARD_METHOD(getCov, coefficientOfVariation, TEST)
  FORWARD_METHOD(getCov, specificVolumeFromCoV, TEST)
  FORWARD_METHOD(getCov, extensionVariance, TEST)
  FORWARD_METHOD(getCov, calculateStDev, TEST)
  FORWARD_METHOD(getCov, evaluateMatInPlace)
  FORWARD_METHOD(getCov, evaluateOneGeneric, TEST)
  FORWARD_METHOD(getCov, evaluateOneIncr, TEST)
  FORWARD_METHOD(getCov, buildVmapOnDbGrid)
  FORWARD_METHOD(getCov, sample)
  FORWARD_METHOD(getCov, sampleUnitary)
  FORWARD_METHOD(getCov, envelop)
  FORWARD_METHOD(getCov, gofToVario, TEST)
  FORWARD_METHOD(getCov, isNoStat)
  FORWARD_METHOD(getCov, manage)
  FORWARD_METHOD(getCov, optimizationPreProcessForData)
  FORWARD_METHOD(getCov, optimizationPostProcess)
  FORWARD_METHOD_NON_CONST(getCov, setOptimEnabled)
  FORWARD_METHOD_NON_CONST(getCov, attachNoStatDb)
  FORWARD_METHOD_NON_CONST(getCov, makeStationary)

  FORWARD_METHOD_NON_CONST(_getCovModify, setContext)

  // Forwarding the methods from _driftList
  
  FORWARD_METHOD(getDriftList, getDrift)
  FORWARD_METHOD(getDriftList, computeDrift, TEST)
  FORWARD_METHOD(getDriftList, evalDriftValue, TEST)
  FORWARD_METHOD(getDriftList, evalDriftMat)
  FORWARD_METHOD(getDriftList, evalDriftMatInPlace)
  FORWARD_METHOD(getDriftList, evalDriftMatByRanks)
  FORWARD_METHOD(getDriftList, evalMeanVecByRanks)
  FORWARD_METHOD(getDriftList, evalDriftMatByRanksInPlace)
  FORWARD_METHOD(getDriftList, evalDriftMatByTargetInPlace)
  FORWARD_METHOD(getDriftList, getNDrift)
  FORWARD_METHOD(getDriftList, getNDriftEquation)
  FORWARD_METHOD(getDriftList, getNExtDrift)
  FORWARD_METHOD(getDriftList, isFlagLinked)
  FORWARD_METHOD(getDriftList, getDriftMaxIRFOrder,-1)
  FORWARD_METHOD(getDriftList, getRankFex)
  FORWARD_METHOD(getDriftList, isDriftSampleDefined)
  FORWARD_METHOD(getDriftList, isDriftFiltered)
  FORWARD_METHOD(getDriftList, isDriftDefined)
  FORWARD_METHOD(getDriftList, isDriftDifferentDefined)
  FORWARD_METHOD(getDriftList, getDrifts)
  FORWARD_METHOD(getDriftList, evalDrift, TEST)
  FORWARD_METHOD(getDriftList, evalDriftBySample)
  FORWARD_METHOD(getDriftList, evalDriftBySampleInPlace)
  FORWARD_METHOD(getDriftList, evalDriftCoef)
  FORWARD_METHOD(getDriftList, hasDrift, false)

  FORWARD_METHOD(getDriftList, getMean, TEST)
  FORWARD_METHOD(getDriftList, getMeans)
  FORWARD_METHOD(getDriftList, evalDriftVarCoef,TEST)
  FORWARD_METHOD(getDriftList, evalDriftVarCoefs)

  FORWARD_METHOD_NON_CONST(_getDriftListModify, setFlagLinked)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, setBetaHat)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, setFiltered)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, delDrift)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, delAllDrifts)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, copyCovContext)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, setMeans)
  FORWARD_METHOD_NON_CONST(_getDriftListModify, setMean)
  
  // Forwarding the methods from _ctxt
  FORWARD_METHOD(getContext, getNVar, -1)
  FORWARD_METHOD(getContext, getNDim, -1)
  FORWARD_METHOD(getContext, getSpace)

  FORWARD_METHOD(getContext, getCovar0)
  FORWARD_METHOD(getContext, getField, TEST)

  FORWARD_METHOD_NON_CONST(_getContextModify, setField)
  FORWARD_METHOD_NON_CONST(_getContextModify, setCovar0s)
  FORWARD_METHOD_NON_CONST(_getContextModify, setCovar0)
  
  void setField(double field);
  bool isValid() const;

  void setCov(ACov* cova);
  
  void setDriftList(const DriftList* driftlist);
  void setDriftIRF(int order = 0, int nfex = 0);
  void addDrift(const ADrift* drift); // TODO: check that the same driftM has not been already defined
  void setDrifts(const VectorString& driftSymbols);

  void initParams();

  #ifndef SWIG
  std::shared_ptr<ListParams> generateListParams() const;
  #endif
  void updateModel();
  double computeLogLikelihood(const Db* db, bool verbose = false);
  double evalGradParam(int iparam, SpacePoint& p1, SpacePoint& p2,int ivar = 0, int jvar = 0);
  void fitNew(const Db* db = nullptr,
              Vario* vario = nullptr,
              const DbGrid* dbmap = nullptr,
              Constraints* constraints = nullptr,
              const ModelOptimParam& mop = ModelOptimParam(),
              int nb_neighVecchia = 30,
              bool verbose = false);

private:
  virtual bool _isValid() const;

protected:               // TODO : pass into private to finish clean
  ACov* _cova;           /* Generic Covariance structure */
  std::vector<std::function<double(double)>> _gradFuncs;
  DriftList* _driftList; /* Series of Drift functions */
  CovContext _ctxt;      /* Context */
};

GSTLEARN_EXPORT int computeCovMatSVCLHSInPlace(MatrixSymmetric& cov,
                                               const MatrixSymmetric& Sigma,
                                               const MatrixDense& F1,
                                               int type = 1,
                                               int idx  = 0);
GSTLEARN_EXPORT int computeCovMatSVCRHSInPlace(MatrixDense& cov,
                                               const MatrixSymmetric& Sigma,
                                               const MatrixDense& F1,
                                               const MatrixDense& F2,
                                               int type1 = 1,
                                               int idx1  = 0,
                                               int type2 = 1,
                                               int idx2  = 0);
GSTLEARN_EXPORT int computeDriftMatSVCRHSInPlace(MatrixDense& mat,
                                                 const MatrixDense& F,
                                                 int type                 = 1,
                                                 int idx                  = 0,
                                                 bool flagCenteredFactors = true);
