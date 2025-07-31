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
#include "API/SPDE.hpp"
#include "API/SPDEParam.hpp"
#include "Basic/Law.hpp"
#include "Basic/NamingConvention.hpp"
#include "Basic/VectorHelper.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovAniso.hpp"
#include "Db/Db.hpp"
#include "Enum/ECov.hpp"
#include "LinearOp/InvNuggetOp.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "LinearOp/PrecisionOpMulti.hpp"
#include "LinearOp/PrecisionOpMultiMatrix.hpp"
#include "LinearOp/ProjMatrix.hpp"
#include "LinearOp/ProjMultiMatrix.hpp"
#include "LinearOp/SPDEOp.hpp"
#include "LinearOp/SPDEOpMatrix.hpp"
#include "LinearOp/ShiftOpMatrix.hpp"
#include "LithoRule/RuleProp.hpp"
#include "Mesh/MeshETurbo.hpp"
#include "Model/Model.hpp"
#include "geoslib_define.h"

#include <cmath>

namespace gstlrn
{

/**
 * The class constructor with the following arguments:
 *
 * @param model  This compulsory argument is a LMC of Matern's (or Markov?) basic structures with possibly a nugget effect
 * @param useCholesky Define the choice regarding Cholesky
 * @param params Set of SPDE parameters
 */
SPDE::SPDE(Model* model,
           int useCholesky,
           const SPDEParam& params)
  : _dbin()
  , _dbout()
  , _model(model)
  , _flagCholesky(false)
  , _driftCoeffs()
  , _createMeshesK(false)
  , _meshesK()
  , _createMeshesS(false)
  , _meshesS()
  , _createAinK(false)
  , _AinK(nullptr)
  , _createAinS(false)
  , _AinS(nullptr)
  , _createAoutK(false)
  , _AoutK(nullptr)
  , _createAoutS(false)
  , _AoutS(nullptr)
  , _QopK(nullptr)
  , _QopS(nullptr)
  , _Qom(nullptr)
  , _invnoiseobj(nullptr)
  , _params(params)
{
  _defineFlagCholesky(useCholesky, model);
}

SPDE::~SPDE()
{
  if (_createMeshesK && !_meshesK.empty())
  {
    for (int i = 0, n = _meshesK.size(); i < n; i++)
      delete _meshesK[i];
  }
  if (_createMeshesS && !_meshesS.empty())
  {
    for (int i = 0, n = _meshesS.size(); i < n; i++)
      delete _meshesS[i];
  }

  if (_createAinK) delete _AinK;
  if (_createAinS) delete _AinS;
  if (_createAoutK) delete _AoutK;
  if (_createAoutS) delete _AoutS;
  delete _Qom;
  delete _QopK;
  delete _QopS;
  delete _invnoiseobj;
}

/**
 * Define if Cholesky must be used or not
 * @param useCholesky: 1 for YES; 0 for No; -1: set optimal default (according to NDim)
 * @param model: The ModelGeneric to use
 * @param verbose: Verbose flag
 */
void SPDE::_defineFlagCholesky(int useCholesky, const Model* model, bool verbose)
{
  if (useCholesky == -1)
  {
    _flagCholesky = model->getNDim() == 2;
  }
  else
    _flagCholesky = (useCholesky == 1);

  if (verbose)
  {
    if (_flagCholesky)
      message("- Choice for the Cholesky option = ON");
    else
      message("- Choice for the Cholesky option = OFF");
    if (useCholesky == -1)
      message(" (Automatic setting)\n");
    else
      message("\n");
  }
}

VectorMeshes SPDE::_defineMeshFromDbs(bool flagKrige)
{
  // Create the domain (by merging dbin and dbout)
  bool isBuilt     = false;
  const Db* domain = Db::coverSeveralDbs(_dbin, _dbout, &isBuilt);

  int refine  = (flagKrige) ? _params.getRefineK() : _params.getRefineS();
  int ncovtot = _model->getNCov(false);
  int ncov    = _model->getNCov(true);
  VectorMeshes meshes(ncov);

  // Create as many meshes as they are structures (nugget excluded)
  for (int jcov = 0, icov = 0; jcov < ncovtot; jcov++)
  {
    if (_model->getCovType(jcov) == ECov::NUGGET) continue;
    const CovAniso* cova = _model->getCovAniso(jcov);
    meshes[icov++]       = MeshETurbo::createFromCova(*cova, domain, refine,
                                                      _params.getBorder(),
                                                      _params.isPolarized(), true,
                                                      _params.getNxMax());
  }
  if (isBuilt) delete domain;

  return meshes;
}

void SPDE::_printMeshesDetails(const VectorMeshes& meshes)
{
  int nmesh = (int)meshes.size();
  for (int imesh = 0; imesh < nmesh; imesh++)
  {
    message("- Mesh #%d/%d: NMeshes = %d, NApices = %d\n",
            imesh + 1, nmesh, meshes[imesh]->getNMeshes(), meshes[imesh]->getNApices());
  }
}

/**
 * @brief Create the set of Meshes needed within SPDE
 *
 * @param flagKrige True if this method is used for Kriging (rather than Simulation)
 * @param meshesIn Input Set of Mesh (one per structure, nugget excluded)
 * @param verbose Verbose flag
 * @return int Error return code
 */
int SPDE::_defineMesh(bool flagKrige, const VectorMeshes& meshesIn, bool verbose)
{
  auto& meshes     = flagKrige ? _meshesK : _meshesS;
  auto& createMesh = flagKrige ? _createMeshesK : _createMeshesS;
  int ncov         = _model->getNCov(true);
  auto nmesh       = meshesIn.size();

  createMesh = false;
  if (nmesh == 0)
  {
    if (!flagKrige && _meshesS.empty() && !_meshesK.empty())
    {
      // Particular case of Simulations: if the corresponding mesh is not defined
      // the meshing dedicated to Kriging is used instead.
      meshes = _meshesK;
      if (verbose) message("Copy Meshings from Kriging meshings\n");
    }
    if (meshes.empty())
    {
      meshes     = _defineMeshFromDbs(flagKrige);
      createMesh = true;
      if (verbose)
      {
        message("Creating Meshes covering the available Db(s)\n");
        _printMeshesDetails(meshes);
      }
    }
  }
  else if (nmesh == 1 && (int)nmesh != ncov)
  {
    meshes = meshesIn;
    // Particular case of a single mesh: simply duplicate it (without creating new contents)
    meshes.resize(ncov);
    for (int icov = 0; icov < ncov; icov++)
      meshes[icov] = meshesIn[0];
    if (verbose)
    {
      message("Duplicating the Input Mesh for each one of the %d covariance(s)\n", ncov);
      _printMeshesDetails(meshes);
    }
  }
  else if ((int)nmesh == ncov)
  {
    // If the dimensions 'nmesh' and 'ncov' match: copy pointer of the input VectorMeshes
    meshes = meshesIn;
  }
  else
  {
    messerr("Argument 'meshes' contains %d items", nmesh);
    messerr("whereas the number of structures (nugget excluded) is %d", ncov);
    return 1;
  }

  return 0;
}

/**
 * @brief Create the set of Projections needed within SPDE
 *
 * @param flagIn True to define the projection system for Input (resp. Output) File
 * @param flagKrige True for Kriging; False for Simulations
 * @param projIn In Projection system (optional)
 * @param verbose Verbose flag
 * @return Pointer to the ProjMultiMatrix used in output (or nullptr)
 */
int SPDE::_defineProjection(bool flagIn,
                            bool flagKrige,
                            const ProjMultiMatrix* projIn,
                            bool verbose)
{
  auto& meshes = flagKrige ? _meshesK : _meshesS;
  const ProjMultiMatrix** projPtr =
    flagIn
      ? (flagKrige ? &_AinK : &_AinS)
      : (flagKrige ? &_AoutK : &_AoutS);
  auto& createProj = flagIn ? (flagKrige ? _createAinK : _createAinS) : (flagKrige ? _createAoutK : _createAoutS);
  const auto* db   = flagIn ? _dbin : _dbout;

  // If no mesh is provided, an empty ProjMultiMatrix is returned

  createProj = false;
  if (meshes.empty()) return 0;

  // Get the number of structures
  int ncov = _model->getNCov(true);
  int nvar = _model->getNVar();

  // Case where the projection matrix has not been provided, create it
  if (projIn == nullptr)
  {
    *projPtr   = ProjMultiMatrix::createFromDbAndMeshes(db, meshes, ncov, nvar, flagIn);
    createProj = true;
    if (verbose)
    {
      if (flagIn)
      {
        if (flagKrige)
          message("Creating Projection from 'dbin' and Mesh(es) for Kriging\n");
        else
          message("Creating Projection from 'dbin' and Mesh(es) for Simulation\n");
      }
      else
      {
        if (flagKrige)
          message("Creating Projection from 'dbout' and Mesh(es) for Kriging\n");
        else
          message("Creating Projection from 'dbout' and Mesh(es) for Simulation\n");
      }
    }
    return 0;
  }

  // The projection matrix is provided: check consistency of its dimensions
  int npoint  = projIn->getNPoint();
  int napices = projIn->getNApex();

  int ncol = 0;
  for (int icov = 0; icov < ncov; icov++)
    for (int ivar = 0; ivar < nvar; ivar++)
      ncol += meshes[icov]->getNApices();
  if (ncol != napices)
  {
    messerr("Number of Columns of 'projm' (%d) is not correct", napices);
    messerr("- Number of covariances = %d", ncov);
    messerr("- Number of variables = %d", nvar);
    for (int icov = 0; icov < nvar; icov++)
      messerr("- Number of apices for Meshing (%d) = %d",
              icov + 1, meshes[icov]->getNApices());
    return 1;
  }

  int nrow = 0;
  for (int ivar = 0; ivar < nvar; ivar++)
    nrow += db->getNSampleActiveAndDefined(ivar);
  if (nrow != npoint)
  {
    messerr("Number of Rows of 'projm' (%d) is not correct", npoint);
    messerr("- Number of variables = %d", nvar);
    for (int ivar = 0; ivar < nvar; ivar++)
      messerr("- Number of samples for variable %d = %d\n",
              ivar, db->getNSampleActiveAndDefined(ivar));
    return 1;
  }

  // All checks are correct: copy the pointer
  *projPtr = projIn;
  if (verbose) message("Copy the input 'projection'\n");

  return 0;
}

SPDEOp* SPDE::defineShiftOperator(bool flagSimu, bool verbose)
{
  SPDEOp* spdeop = nullptr;

  _invnoiseobj = new InvNuggetOp(_dbin, _model, _params, !_flagCholesky);

  if (_flagCholesky)
  {
    if (!_meshesK.empty())
    {
      _Qom   = new PrecisionOpMultiMatrix(_model, _meshesK);
      spdeop = new SPDEOpMatrix(_Qom, _AinK, _invnoiseobj, _AoutK);
    }
    else
    {
      _Qom   = new PrecisionOpMultiMatrix(_model, _meshesS);
      spdeop = new SPDEOpMatrix(_Qom, _AinK, _invnoiseobj, _AoutS);
    }
  }
  else
  {
    _QopK = new PrecisionOpMulti(_model, _meshesK, _params.getUseStencil());
    if (flagSimu)
      _QopS = new PrecisionOpMulti(_model, _meshesS, _params.getUseStencil());

    spdeop = new SPDEOp(_QopK, _AinK, _invnoiseobj, _QopS, _AinS, _AoutK, _AoutS);
    spdeop->setMaxIterations(_params.getCGparams().getNIterMax());
    spdeop->setTolerance(_params.getCGparams().getEps());
  }

  spdeop->setVerbose(verbose);
  return spdeop;
}

int SPDE::defineMeshes(bool flagSimu,
                       const VectorMeshes& meshesK,
                       const VectorMeshes& meshesS,
                       bool verbose)
{
  if (_dbin != nullptr)
  {
    if (_defineMesh(true, meshesK, verbose)) return 1;
  }
  if (flagSimu)
  {
    if (_defineMesh(false, meshesS, verbose)) return 1;
  }
  return 0;
}

int SPDE::defineProjections(bool flagSimu,
                            bool flagCond,
                            const ProjMultiMatrix* projInK,
                            const ProjMultiMatrix* projInS,
                            bool verbose)
{
  if (flagCond && _dbin != nullptr)
  {
    // Projection of Data for Kriging
    if (_defineProjection(true, true, projInK, verbose)) return 1;
    if (flagSimu)
    {
      // Projection of Data for Simulation
      if (_defineProjection(true, false, projInS, verbose)) return 1;
    }
  }

  if (flagCond && _dbout != nullptr)
  {
    // Projection of Target for Kriging
    if (_defineProjection(false, true, nullptr, verbose)) return 1;
  }
  if (flagSimu)
  {
    // Projection of Target for Simulation
    if (_defineProjection(false, false, nullptr, verbose)) return 1;
  }
  return 0;
}

int SPDE::centerDataByDriftInPlace(const SPDEOp* spdeop, VectorDouble& Z, bool verbose)
{
  if (Z.empty()) return 1;
  _driftCoeffs.clear();
  bool mustEvaluateDrift = _model->getNDrift() > 0;

  if (mustEvaluateDrift)
  {
    MatrixDense driftMat = _model->evalDriftMatByRanks(_dbin);
    _driftCoeffs         = spdeop->computeDriftCoeffs(Z, driftMat);
    ASPDEOp::centerDataByDriftMat(Z, driftMat, _driftCoeffs);
    if (verbose) VH::dump("Drift coefficients = ", _driftCoeffs);
  }
  else
  {
    VectorDouble meanVec = _model->evalMeanVecByRanks(_dbin);
    ASPDEOp::centerDataByMeanVec(Z, meanVec);
  }
  return 0;
}

void SPDE::uncenterResultByDriftInPlace(VectorDouble& result)
{
  int nvar               = _model->getNVar();
  int nech               = _dbout->getNSample();
  int nechred            = _dbout->getNSample(true);
  bool mustEvaluateDrift = (_model->getNDrift() > 0);

  // Loop on the samples of the output file

  double value;
  int ncols = 0;
  MatrixDense mat;
  for (int jech = 0, iech = 0; jech < nech; jech++)
  {
    if (!_dbout->isActive(jech)) continue;

    if (mustEvaluateDrift)
    {
      if (_driftCoeffs.empty()) continue;
      (void)_model->evalDriftMatByTargetInPlace(mat, _dbout, jech);
      ncols = mat.getNCols();
    }

    // Loop on the variables

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      value = 0.;
      if (mustEvaluateDrift)
      {
        for (int icol = 0; icol < ncols; icol++)
          value += mat.getValue(ivar, icol) * _driftCoeffs[icol];
      }
      else
        value = _model->getMean(ivar);
      result[ivar * nechred + iech] += value;
    }
    iech++;
  }
}

/**
 * @brief Simulate nugget component and add it to one multivariate simulation
 *
 * @param result The array containing ONE multivariate simulation
 *
 * @note The nugget component is added to each variable with its correct variance
 * @note but not accounting for the dependency across different variables.
 */
void SPDE::addNuggetToResult(VectorDouble& result)
{
  int rankNugget = _model->getRankNugget();
  if (rankNugget < 0) return;
  int nvar = _model->getNVar();
  for (int ivar = 0, ecr = 0; ivar < nvar; ivar++)
  {
    double nugget = _model->getSill(rankNugget, ivar, ivar);
    for (int iech = 0, nech = (int)result.size(); iech < nech; iech++, ecr++)
      result[iech] += law_gaussian(0., sqrt(nugget));
  }
}

/**
 * Derive the global trend in the SPDE framework
 *
 * @param dbin Input Db (must contain the variable to be estimated)
 * @param dbout Output Db where the estimation must be performed
 * @param model Model definition
 * @param useCholesky Define the choice regarding Cholesky (see _defineCholesky)
 * @param meshesK Meshes description (optional)
 * @param projInK Matrix of projection (optional)
 * @param params Set of SPDE parameters
 * @param verbose Verbose flag
 *
 * @return Returned vector
 */
VectorDouble trendSPDE(Db* dbin,
                       Db* dbout,
                       Model* model,
                       int useCholesky,
                       const VectorMeshes& meshesK,
                       const ProjMultiMatrix* projInK,
                       const SPDEParam& params,
                       bool verbose)
{
  // Preliminary checks
  if (dbin == nullptr) return 1;
  if (dbout == nullptr) return 1;
  if (model == nullptr) return 1;

  // Instantiate SPDE class
  SPDE spde(model, useCholesky, params);
  spde.setDbin(dbin);
  spde.setDbout(dbout);
  if (verbose) mestitle(1, "Trend in SPDE framework (Cholesky=%d)", (int)spde.getFlagCholesky());

  // Define Meshes
  if (spde.defineMeshes(false, meshesK, VectorMeshes(), verbose)) return 1;

  // Define projections
  if (spde.defineProjections(false, true, projInK, nullptr, verbose)) return 1;

  // Define the Shift operator
  auto* spdeop = spde.defineShiftOperator(false);

  // Read information from the input Db and center it
  VectorDouble Z = dbin->getColumnsActiveAndDefined(ELoc::Z);
  if (spde.centerDataByDriftInPlace(spdeop, Z, verbose)) return 1;

  delete spdeop;
  return spde.getDriftCoefficients();
}

/**
 * Perform the estimation by KRIGING under the SPDE framework
 *
 * @param dbin Input Db (must contain the variable to be estimated)
 * @param dbout Output Db where the estimation must be performed
 * @param model Model definition
 * @param flag_est True for the estimation
 * @param flag_std True for the standard deviation of estimation error
 * @param useCholesky Define the choice regarding Cholesky (see _defineCholesky)
 * @param meshesK Meshes description (optional)
 * @param projInK Matrix of projection (optional)
 * @param meshesS Meshes used for Variance calulcation (optional)
 * @param projInS Matrix of projection used for Variance calculation (optional)
 * @param params Set of SPDE parameters
 * @param verbose Verbose flag
 * @param namconv Naming convention
 *
 * @return Returned vector
 *
 * @remark Algorithm for 'meshesK' or 'mesheS' and 'projInK' or 'projInS':
 * - Each one of the previous arguments is considered individually and sequentially.
 * - For 'meshes' in general ('meshesK' or 'meshesS'):
 *   - If it is not defined, it is created from 'model' and so as to cover both 'dbin' and 'dbout' (if provided).
 *   - If is already exist, the number of meshes must be equal:
 *     - either to 1: the same mesh is used for all structures
 *     - or to the number of structures (nugget discarded)
 *   - Otherwise an error is raised
 * - For 'projIn' in general ('projInK' or 'projInS'):
 *   - If it is not defined, it is created from 'meshesK'
 *   - If it is already exist, the number of projectors must be equal to the number of meshes
 *   - Otherwise an error is raised
 * - If 'mesheS' does not exist, 'meshesK' is used instead
 * - If 'projInS' does not exist, 'projInK' is used instead
 *
 * 'meshesS' and 'projInS' are used only if 'flag_std' is true and useCholesky=0
 */
int krigingSPDE(Db* dbin,
                Db* dbout,
                Model* model,
                bool flag_est,
                bool flag_std,
                int useCholesky,
                const VectorMeshes& meshesK,
                const ProjMultiMatrix* projInK,
                const VectorMeshes& meshesS,
                const ProjMultiMatrix* projInS,
                const SPDEParam& params,
                bool verbose,
                const NamingConvention& namconv)
{
  // Preliminary checks
  if (dbin == nullptr) return 1;
  if (dbout == nullptr) return 1;
  if (model == nullptr) return 1;

  // Instantiate SPDE class
  SPDE spde(model, useCholesky, params);
  bool flagSimu = flag_std && !spde.getFlagCholesky();
  spde.setDbin(dbin);
  spde.setDbout(dbout);
  if (verbose) mestitle(1, "Kriging in SPDE framework (Cholesky=%d)", (int)spde.getFlagCholesky());

  // Define Meshes
  if (spde.defineMeshes(flagSimu, meshesK, meshesS, verbose)) return 1;

  // Define projections
  if (spde.defineProjections(flagSimu, true, projInK, projInS, verbose)) return 1;

  // Define the Shift operator
  auto* spdeop = spde.defineShiftOperator(flagSimu);

  // Read information from the input Db and center it
  VectorDouble Z = dbin->getColumnsActiveAndDefined(ELoc::Z);
  if (spde.centerDataByDriftInPlace(spdeop, Z, verbose)) return 1;

  // Performing the task and storing results in 'dbout'
  // This is performed in ONE step to avoid additional core allocation
  int nvar = model->getNVar();
  VectorDouble result;
  if (flag_est)
  {
    result = spdeop->kriging(Z);
    spde.uncenterResultByDriftInPlace(result);
    int iuid = dbout->addColumns(result, "estim", ELoc::Z, 0, true, 0., nvar);
    namconv.setNamesAndLocators(dbin, VectorString(), ELoc::Z, nvar, dbout, iuid,
                                "estim");
  }
  if (flag_std)
  {
    result   = spdeop->stdev(Z, spde.getNMC(), spde.getSeed());
    int iuid = dbout->addColumns(result, "stdev", ELoc::UNKNOWN, 0, true, 0., nvar);
    namconv.setNamesAndLocators(dbin, VectorString(), ELoc::Z, nvar, dbout, iuid,
                                "stdev");
  }

  delete spdeop;
  return 0;
}

/**
 * Perform the conditional SIMULATIONs in the SPDE framework
 *
 * @param dbin Input Db (variable to be estimated). Only for conditional simulations
 * @param dbout Output Db where the estimation must be performed
 * @param model Model definition
 * @param nbsimu Number of simulations
 * @param useCholesky Define the choice regarding Cholesky (see _defineCholesky)
 * @param meshesK Meshes used for Kriging (optional)
 * @param projInK Matrix of projection used for Kriging (optional)
 * @param meshesS Meshes used for Simulations (optional)
 * @param projInS Matrix of projection used for Simulations (optional)
 * @param params Set of SPDE parameters
 * @param verbose Verbose flag
 * @param namconv see NamingConvention
 *
 * @return Error returned code
 *
 * @remark Algorithm for 'meshesK', 'projInK', 'meshesS' and 'projInS':
 * - Each one of the previous arguments is considered individually and sequentially.
 * - For 'meshes' in general ('meshesK' or 'meshesS'):
 *   - If it is not defined, it is created from 'model' and so as to cover both 'dbin' and 'dbout' (if provided).
 *   - If is already exist, the number of meshes must be equal:
 *     - either to 1: the same mesh is used for all structures
 *     - or to the number of structures (nugget discarded)
 *   - Otherwise an error is raised
 * - For 'projIn' in general ('projInK' or 'projInS'):
 *   - If it is not defined, it is created from 'meshes'
 *   - If it is already exist, the number of projectors must be equal to the number of meshes
 *   - Otherwise an error is raised
 * - If 'mesheS' does not exist, 'meshesK' is used instead
 * - If 'projInS' does not exist, 'projInK' is used instead
 */
int simulateSPDE(Db* dbin,
                 Db* dbout,
                 Model* model,
                 int nbsimu,
                 int useCholesky,
                 const VectorMeshes& meshesK,
                 const ProjMultiMatrix* projInK,
                 const VectorMeshes& meshesS,
                 const ProjMultiMatrix* projInS,
                 const SPDEParam& params,
                 bool verbose,
                 const NamingConvention& namconv)
{
  if (dbout == nullptr) return 1;
  if (model == nullptr) return 1;
  bool flagCond = (dbin != nullptr);

  // Instantiate SPDE class
  SPDE spde(model, useCholesky, params);
  spde.setDbin(dbin);
  spde.setDbout(dbout);
  if (verbose) mestitle(1, "Simulation in SPDE framework (cond=%d, Cholesky=%d)",
                        (int)flagCond, (int)spde.getFlagCholesky());

  // Define Meshes
  if (spde.defineMeshes(true, meshesK, meshesS, verbose)) return 1;

  // Define
  if (spde.defineProjections(true, flagCond, projInK, projInS, verbose)) return 1;

  // Define the Shift operator
  auto* spdeop = spde.defineShiftOperator(true);

  VectorDouble Z;
  VectorDouble driftCoeffs;
  if (flagCond)
  {
    // Read information from the input Db and center it
    Z = dbin->getColumnsActiveAndDefined(ELoc::Z);
    if (spde.centerDataByDriftInPlace(spdeop, Z, verbose)) return 1;
  }

  // Perform the Simulation and storage.
  // All is done in ONE step to avoid additional storage
  int nvar    = model->getNVar();
  int iuid    = dbout->addColumnsByConstant(nvar * nbsimu);
  int nechred = dbout->getNSample(true);
  VectorDouble local(nechred);
  VectorDouble result;

  for (int isimu = 0; isimu < nbsimu; isimu++)
  {
    result = (flagCond) ? spdeop->simCond(Z) : spdeop->simNonCond();
    spde.addNuggetToResult(result);
    spde.uncenterResultByDriftInPlace(result);

    for (int ivar = 0; ivar < nvar; ivar++)
    {
      VH::extractInPlace(result, local, ivar * nechred);
      int juid = iuid + ivar * nbsimu + isimu;
      dbout->setColumnByUID(local, juid, true);
    }
  }
  namconv.setNamesAndLocators(dbin, VectorString(), ELoc::Z, nvar, dbout, iuid,
                              "", nbsimu);

  delete spdeop;
  return 0;
}

/**
 * Perform the conditional SIMULATIONs using PGS technique in the SPDE framework
 *
 * @param dbin Input Db (variable to be estimated). Only for conditional simulations
 * @param dbout Output Db where the estimation must be performed
 * @param model Model definition
 * @param ruleprop RuleProp structure describing the Rule and the Proportions
 * @param nbsimu Number of simulations
 * @param useCholesky Define the choice regarding Cholesky (see _defineCholesky)
 * @param meshesK Meshes used for Kriging (optional)
 * @param projInK Matrix of projection used for Kriging (optional)
 * @param meshesS Meshes used for Simulations (optional)
 * @param projInS Matrix of projection used for Simulations (optional)
 * @param params Set of SPDE parameters
 * @param verbose Verbose flag
 * @param namconv see NamingConvention
 *
 * @return Error returned code
 *
 * @remark Algorithm for 'meshesK', 'projInK', 'meshesS' and 'projInS':
 * - Each one of the previous arguments is considered individually and sequentially.
 * - For 'meshes' in general ('meshesK' or 'meshesS'):
 *   - If it is not defined, it is created from 'model' and so as to cover both 'dbin' and 'dbout' (if provided).
 *   - If is already exist, the number of meshes must be equal:
 *     - either to 1: the same mesh is used for all structures
 *     - or to the number of structures (nugget discarded)
 *   - Otherwise an error is raised
 * - For 'projIn' in general ('projInK' or 'projInS'):
 *   - If it is not defined, it is created from 'meshes'
 *   - If it is already exist, the number of projectors must be equal to the number of meshes
 *   - Otherwise an error is raised
 * - If 'mesheS' does not exist, 'meshesK' is used instead
 * - If 'projInS' does not exist, 'projInK' is used instead
 */
int simPGSSPDE(Db* dbin,
               Db* dbout,
               Model* model,
               const RuleProp& ruleprop,
               int nbsimu,
               int useCholesky,
               const VectorMeshes& meshesK,
               const ProjMultiMatrix* projInK,
               const VectorMeshes& meshesS,
               const ProjMultiMatrix* projInS,
               const SPDEParam& params,
               bool verbose,
               const NamingConvention& namconv)
{
  if (dbout == nullptr) return 1;
  if (model == nullptr) return 1;
  bool flagCond = (dbin != nullptr);

  // Instantiate SPDE class
  SPDE spde(model, useCholesky, params);
  spde.setDbin(dbin);
  spde.setDbout(dbout);
  if (verbose) mestitle(1, "PluriGaussian Simulation in SPDE framework (cond=%d, Cholesky=%d)",
                        (int)flagCond, (int)spde.getFlagCholesky());

  // Define Meshes
  if (spde.defineMeshes(true, meshesK, meshesS, verbose)) return 1;

  // Define projections
  if (spde.defineProjections(true, flagCond, projInK, projInS, verbose)) return 1;

  // Define the Shift operator
  auto* spdeop = spde.defineShiftOperator(true);

  VectorDouble Z;
  if (flagCond)
  {
    // Read information from the input Db and center it
    Z = dbin->getColumnsActiveAndDefined(ELoc::Z);
    if (spde.centerDataByDriftInPlace(spdeop, Z, verbose)) return 1;
  }

  // Perform the Simulation and storage.
  // All is done in ONE step to avoid additional storage
  int nvar = model->getNVar();

  int nechred = dbout->getNSample(true);
  VectorDouble local(nechred);
  VectorDouble result;

  if (flagCond)
    ruleprop.categoryToThresh(dbin);

  // Loop on the simulations
  for (int isimu = 0; isimu < nbsimu; isimu++)
  {
    int iuid = dbout->addColumnsByConstant(nvar);

    result = (flagCond) ? spdeop->simCond(Z) : spdeop->simNonCond();
    spde.addNuggetToResult(result);
    spde.uncenterResultByDriftInPlace(result);

    // Loop on the variables
    for (int ivar = 0; ivar < nvar; ivar++)
    {
      VH::extractInPlace(result, local, ivar * nechred);
      dbout->setColumnByUID(local, iuid + ivar, true);
      dbout->setLocatorByUID(iuid + ivar, ELoc::SIMU, ivar);
    }

    // Convert the resulting simulation into categories
    ruleprop.gaussToCategory(dbout, namconv);
    dbout->deleteColumnsByUID(VH::sequence(nvar, iuid));
  }

  delete spdeop;
  return 0;
}

/**
 * Calculate the Log-Likelihood under the SPDE framework
 *
 * @param dbin Input Db (must contain the variable to be estimated)
 * @param model Model definition
 * @param useCholesky Define the choice regarding Cholesky (see _defineCholesky)
 * @param meshes Meshes description (optional)
 * @param projIn Matrix of projection (optional)
 * @param params Set of SPDE parameters
 * @param verbose True for verbose output
 * @return Returned value
 */
double logLikelihoodSPDE(Db* dbin,
                         Model* model,
                         int useCholesky,
                         const VectorMeshes& meshes,
                         const ProjMultiMatrix* projIn,
                         const SPDEParam& params,
                         bool verbose)
{
  if (dbin == nullptr) return 1;
  if (model == nullptr) return 1;
  if (dbin->getNLoc(ELoc::Z) != 1)
  {
    messerr("'dbin' must contain ONE variable (Z locator)");
    return 1;
  }

  // Instantiate SPDE class
  SPDE spde(model, useCholesky, params);
  spde.setDbin(dbin);
  if (verbose) mestitle(1, "Log-likelhood calculation in SPDE framework( Cholesky=%d)",
                        (int)spde.getFlagCholesky());

  // Define Meshes
  if (spde.defineMeshes(true, meshes, VectorMeshes(), verbose)) return 1;

  // Define projections
  if (spde.defineProjections(false, true, projIn, nullptr, verbose)) return 1;

  // Define the Shift operator
  auto* spdeop = spde.defineShiftOperator(false, verbose);

  // Read information from the input Db and center it
  VectorDouble Z = dbin->getColumnsActiveAndDefined(ELoc::Z);
  if (spde.centerDataByDriftInPlace(spdeop, Z, verbose)) return 1;

  // Performing the task
  int size       = (int)Z.size();
  double logdet  = spdeop->computeTotalLogDet(spde.getNMC(), spde.getSeed());
  double quad    = spdeop->computeQuadratic(Z);
  double loglike = TEST;
  if (!FFFF(logdet) && !FFFF(quad))
    loglike = -0.5 * (logdet + quad + size * log(2. * GV_PI));

  if (verbose)
  {
    message("Likelihood calculation:\n");
    message("Nb. active samples = %d\n", size);
    message("Nb. Monte-Carlo    = %d\n", spde.getNMC());
    message("Cholesky           = %d\n", spde.getFlagCholesky());
    message("Log-Determinant    = %lf\n", logdet);
    message("Quadratic term     = %lf\n", quad);
    message("Log-likelihood     = %lf\n", loglike);
  }

  // Cleaning phase
  delete spdeop;
  return loglike;
}

} // namespace gstlrn