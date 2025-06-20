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

#include "Estimation/AModelOptimNew.hpp"
#include "gstlearn_export.hpp"

#include "Basic/VectorNumT.hpp"
#include "Model/ModelOptimParam.hpp"

class Model;
class DbGrid;
class Constraints;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * various parameters of a Covariance part of a Model
 */
class GSTLEARN_EXPORT ModelOptimVMap: public AModelOptimNew
{
public:
  ModelOptimVMap(ModelGeneric* model,
                 Constraints* constraints   = nullptr,
                 const ModelOptimParam& mop = ModelOptimParam());
  ModelOptimVMap(const ModelOptimVMap& m);
  ModelOptimVMap& operator=(const ModelOptimVMap& m);
  virtual ~ModelOptimVMap();

  double computeCost(bool verbose = false) override;

  static ModelOptimVMap* createForOptim(ModelGeneric* model,
                                        const DbGrid* dbmap,
                                        Constraints* constraints   = nullptr,
                                        const ModelOptimParam& mop = ModelOptimParam());

private:
  bool _checkConsistency();
  int  _getDimensions();
  void _allocateInternalArrays();

protected:
  // Model fitting options
  ModelOptimParam _mop;

  // Set of constraints
  Constraints* _constraints;

  // Calculation option
  CovCalcMode _calcmode;

  // Part relative to the Experimental VMap
  const DbGrid* _dbmap;

  // Following members are simply there to accelerate the computation
  VectorInt _indg1;
  VectorInt _indg2;
  int _ndim;
  int _nvar;
  int _nech;
  int _npadir;
};
