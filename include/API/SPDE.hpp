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

#include "API/SPDEParam.hpp"
#include "Basic/NamingConvention.hpp"
#include "LinearOp/InvNuggetOp.hpp"
#include "LinearOp/PrecisionOpMatrix.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"
#include "LinearOp/PrecisionOpMultiMatrix.hpp"
#include "LinearOp/ProjMultiMatrix.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
class ShiftOpMatrix;
class PrecisionOp;
class PrecisionOpMatrix;
class Db;
class DbGrid;
class MeshETurbo;
class Model;
class RuleProp;
class SPDEOp;

/**
 * The SPDE class provides the SPDE implementation of a univariate model defined by
 * the sum of a nugget effect and Matern's models. Its main objectives are:
 * - point kriging with or without linear drifts (SK, OK, KED, UK)
 * - point simulations, conditional or non conditional
 * - evaluation of the log likelihood of the model, in order to estimate the parameter using maximum likelihood
 */
class GSTLEARN_EXPORT SPDE
{
public:
  SPDE(Model* model,
       int useCholesky         = -1,
       const SPDEParam& params = SPDEParam());
  SPDE(const SPDE& r)            = delete;
  SPDE& operator=(const SPDE& r) = delete;
  virtual ~SPDE();

public:
  bool getFlagCholesky() const { return _flagCholesky; }
  void setDbin(const Db* dbin) { _dbin = dbin; };
  void setDbout(const Db* dbout) { _dbout = dbout; };

  const VectorMeshes& getMeshesK() const { return _meshesK; }
  const VectorMeshes& getMeshesS() const { return _meshesS; }
  const ProjMultiMatrix* getAinK() const { return _AinK; }
  const ProjMultiMatrix* getAinS() const { return _AinS; }
  const ProjMultiMatrix* getAoutK() const { return _AoutK; }
  const ProjMultiMatrix* getAoutS() const { return _AoutS; }
  VectorDouble getDriftCoefficients() const { return _driftCoeffs; }

  int getSeed() const { return _params.getSeedMC(); }
  int getNMC() const { return _params.getNMC(); }

  int defineMeshes(bool flagSimu,
                   const VectorMeshes& meshesK,
                   const VectorMeshes& meshesS = VectorMeshes(),
                   bool verbose                = false);
  int defineProjections(bool flagSimu,
                        bool flagCond,
                        const ProjMultiMatrix* projInK,
                        const ProjMultiMatrix* projInS = nullptr,
                        bool verbose                   = false);
  SPDEOp* defineShiftOperator(bool flagSimu, bool verbose = false);
  int centerDataByDriftInPlace(const SPDEOp* spdeop, VectorDouble& Z, bool verbose = false);
  void uncenterResultByDriftInPlace(VectorDouble& result);
  void addNuggetToResult(VectorDouble& result);

private:
  void _defineFlagCholesky(int useCholesky, const Model* model, bool verbose = false);
  VectorMeshes _defineMeshFromDbs(bool flagKrige);
  int _defineMesh(bool flagKrige, const VectorMeshes& meshesIn, bool verbose = false);
  int _defineProjection(bool flagIn,
                        bool flagKrige,
                        const ProjMultiMatrix* projIn,
                        bool verbose = false);
  static void _printMeshesDetails(const VectorMeshes& meshes);

private:
  const Db* _dbin;  // External Pointer
  const Db* _dbout; // External pointer
  Model* _model;    // External pointer
  bool _flagCholesky;
  VectorDouble _driftCoeffs;
  bool _createMeshesK;
  VectorMeshes _meshesK;
  bool _createMeshesS;
  VectorMeshes _meshesS;
  bool _createAinK;
  const ProjMultiMatrix* _AinK;
  bool _createAinS;
  const ProjMultiMatrix* _AinS;
  bool _createAoutK;
  const ProjMultiMatrix* _AoutK;
  bool _createAoutS;
  const ProjMultiMatrix* _AoutS;
  PrecisionOpMulti* _QopK;
  PrecisionOpMulti* _QopS;
  PrecisionOpMultiMatrix* _Qom;
  InvNuggetOp* _invnoiseobj;
  SPDEParam _params;
};

GSTLEARN_EXPORT VectorDouble trendSPDE(Db* dbin,
                                       Db* dbout,
                                       Model* model,
                                       int useCholesky                = -1,
                                       const VectorMeshes& meshesK    = VectorMeshes(),
                                       const ProjMultiMatrix* projInK = nullptr,
                                       const SPDEParam& params        = SPDEParam(),
                                       bool verbose                   = false);
GSTLEARN_EXPORT int krigingSPDE(Db* dbin,
                                Db* dbout,
                                Model* model,
                                bool flag_est                   = true,
                                bool flag_std                   = false,
                                int useCholesky                 = -1,
                                const VectorMeshes& meshesK     = VectorMeshes(),
                                const ProjMultiMatrix* projInK  = nullptr,
                                const VectorMeshes& meshesS     = VectorMeshes(),
                                const ProjMultiMatrix* projInS  = nullptr,
                                const SPDEParam& params         = SPDEParam(),
                                bool verbose                    = false,
                                const NamingConvention& namconv = NamingConvention("KrigingSPDE"));
GSTLEARN_EXPORT int simulateSPDE(Db* dbin,
                                 Db* dbout,
                                 Model* model,
                                 int nbsimu                      = 1,
                                 int useCholesky                 = -1,
                                 const VectorMeshes& meshesK     = VectorMeshes(),
                                 const ProjMultiMatrix* projInK  = nullptr,
                                 const VectorMeshes& meshesS     = VectorMeshes(),
                                 const ProjMultiMatrix* projInS  = nullptr,
                                 const SPDEParam& params         = SPDEParam(),
                                 bool verbose                    = false,
                                 const NamingConvention& namconv = NamingConvention("SimuSPDE"));
GSTLEARN_EXPORT int simPGSSPDE(Db* dbin,
                               Db* dbout,
                               Model* model,
                               const RuleProp& ruleprop,
                               int nbsimu                      = 1,
                               int useCholesky                 = -1,
                               const VectorMeshes& meshesK     = VectorMeshes(),
                               const ProjMultiMatrix* projInK  = nullptr,
                               const VectorMeshes& meshesS     = VectorMeshes(),
                               const ProjMultiMatrix* projInS  = nullptr,
                               const SPDEParam& params         = SPDEParam(),
                               bool verbose                    = false,
                               const NamingConvention& namconv = NamingConvention("SimPGSSPDE"));
GSTLEARN_EXPORT double logLikelihoodSPDE(Db* dbin,
                                         Model* model,
                                         int useCholesky               = -1,
                                         const VectorMeshes& meshes    = VectorMeshes(),
                                         const ProjMultiMatrix* projIn = nullptr,
                                         const SPDEParam& params       = SPDEParam(),
                                         bool verbose                  = false);

} // namespace gstlrn