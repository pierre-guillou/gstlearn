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

#include "Model/AModelOptim.hpp"

class Model;
class Db;

/**
 * \brief
 * Class which, starting from an experimental variogram, enables fitting the
 * various parameters of a Covariance part of a Model
 */
class GSTLEARN_EXPORT ModelOptimLikelihood: public AModelOptim
{
public:
  ModelOptimLikelihood(Model* model);
  ModelOptimLikelihood(const ModelOptimLikelihood& m);
  ModelOptimLikelihood& operator=(const ModelOptimLikelihood& m);
  virtual ~ModelOptimLikelihood();

  int fit(Db* db, bool flagSPDE = false, bool verbose = false);
  int loadEnvironment(Db* db, bool flagSPDE = false, bool verbose = false);

#ifndef SWIG
  static double evalCost(unsigned int nparams,
                         const double* current,
                         double* grad,
                         void* my_func_data);
#endif

private:
  struct Db_Part
  {
    // Use SPDE approach (TRUE); use the Covariance Matrix approach (FALSE°
    // Use SPDE approach (TRUE); use the Covariance Matrix approach (FALSE°
    bool _flagSPDE;
    Db* _db;
  };

  struct AlgorithmLikelihood
  {
    Model_Part& _modelPart;
    Db_Part& _dbPart;
  };
private:
private:
  void _copyDbPart(const Db_Part& dbPart);
  bool _checkConsistency();

private:
  // Part relative to the Db
  Db_Part _dbPart;
};
