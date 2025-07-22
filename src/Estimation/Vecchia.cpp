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
#include "Estimation/Vecchia.hpp"

#include "Basic/VectorHelper.hpp"
#include "Db/Db.hpp"
#include "Estimation/ALikelihood.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "LinearOp/CholeskySparse.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixT.hpp"
#include "Model/ModelGeneric.hpp"
#include "Stats/Classical.hpp"
#include "Tree/Ball.hpp"
#include "geoslib_define.h"

namespace gstlrn
{
Vecchia::Vecchia(ModelGeneric* model,
                 int nb_neigh,
                 const Db* db1,
                 const Db* db2,
                 bool reml)
  : ALikelihood(model, db1, reml)
  , _nbNeigh(nb_neigh)
  , _db1(db1)
  , _db2(db2)
  , _Ranks()
  , _matCov(std::make_shared<MatrixSymmetric>(0))
  , _vectCov()
  , _work()
  , _Y()
  , _DFull()
  , _LFull()
  , _Ndb1(0)
  , _Ntot1(0)
  , _Ntot2(0)
  , _cumulRanks1()
  , _cumulRanks2()
  , _varRanks1()
  , _varRanks2()
{
  setAuthorizedAnalyticalGradients(false);
  _chol = new CholeskyDense();
}

Vecchia::Vecchia(const Vecchia& r)
  : ALikelihood(r)
  , _nbNeigh(r._nbNeigh)
  , _db1(r._db1)
  , _db2(r._db2)
  , _Ranks(r._Ranks)
  , _matCov(r._matCov)
  , _vectCov(r._vectCov)
  , _work(r._work)
  , _Y(r._Y)
  , _DFull(r._DFull)
  , _LFull(r._LFull)
  , _chol(nullptr) // Will be initialized in the constructor
  , _Ndb1(r._Ndb1)
  , _Ntot1(r._Ntot1)
  , _Ntot2(r._Ntot2)
  , _cumulRanks1(r._cumulRanks1)
  , _cumulRanks2(r._cumulRanks2)
  , _varRanks1(r._varRanks1)
  , _varRanks2(r._varRanks2)
{
  _chol = new CholeskyDense(*r._chol);
}
Vecchia& Vecchia::operator=(const Vecchia& r)
{
  if (this != &r)
  {
    ALikelihood::operator=(r);
    _nbNeigh     = r._nbNeigh;
    _db1         = r._db1;
    _db2         = r._db2;
    _Ranks       = r._Ranks;
    _matCov      = r._matCov;
    _vectCov     = r._vectCov;
    _work        = r._work;
    _Y           = r._Y;
    _DFull       = r._DFull;
    _LFull       = r._LFull;
    _Ndb1        = r._Ndb1;
    _Ntot1       = r._Ntot1;
    _Ntot2       = r._Ntot2;
    _cumulRanks1 = r._cumulRanks1;
    _cumulRanks2 = r._cumulRanks2;
    _varRanks1   = r._varRanks1;
    _varRanks2   = r._varRanks2;

    delete _chol;
    _chol = new CholeskyDense(*r._chol);
  }
  return *this;
}

Vecchia::~Vecchia()
{
  delete _chol;
}

void Vecchia::_init(bool verbose)
{
  _Ranks = findNN(_db, nullptr, _nbNeigh + 1, false, verbose);
}

int Vecchia::_getAddressAbsolute(int ip) const
{
  if (ip < _Ndb1)
    return ip;
  return ip - _Ndb1;
}

int Vecchia::_getSampleCase(int ip) const
{
  return (ip < _Ndb1) ? 1 : 2;
}

int Vecchia::_getCase() const
{
  int icase = 0;
  if (_db1 != nullptr) icase = 1;
  if (_db2 != nullptr) icase = 2;
  if (_db1 != nullptr && _db2 != nullptr) icase = 2;
  return icase;
}

int Vecchia::_getAddressInMatrix(int ip, int ivar) const
{
  if (ip < _Ndb1)
    return _cumulRanks1[ivar] + _varRanks1[ivar][ip];
  return _Ntot1 + _cumulRanks2[ivar] + _varRanks2[ivar][ip - _Ndb1];
}

int Vecchia::_buildNeighborhood(const MatrixT<int>& Ranks,
                                int isample,
                                int ivar,
                                int nb_neigh,
                                std::vector<std::array<int, 4>>& neighDescr) const
{
  // Loop on the ranks of the neighboring samples
  int nitems = 0;
  for (int jp = 0; jp < nb_neigh; jp++)
  {
    int ip = Ranks(isample, jp + 1);

    // Discard the missing neighborhood index
    if (IFFFF(ip)) continue;

    // Loop on the variables
    for (int jvar = 0; jvar <= ivar; jvar++)
    {
      neighDescr[nitems][0] = _getSampleCase(ip);            // 1 for Db1; 2 for Db2
      neighDescr[nitems][1] = jvar;                          // Rank of the variable
      neighDescr[nitems][2] = _getAddressInMatrix(ip, jvar); // Position in the matrix
      neighDescr[nitems][3] = _getAddressAbsolute(ip);       // Rank of the sample
      nitems++;
    }
  }

  // For high variable rank, ensure that the same sample for previous variables is selected

  if (ivar > 0)
  {
    int ip0 = Ranks(isample, 0);
    for (int jvar = 0; jvar < ivar; jvar++)
    {
      neighDescr[nitems][0] = _getSampleCase(ip0);            // 1 for Db1; 2 for Db2
      neighDescr[nitems][1] = jvar;                           // Rank of the variable
      neighDescr[nitems][2] = _getAddressInMatrix(ip0, jvar); // Position in the matrix
      neighDescr[nitems][3] = _getAddressAbsolute(ip0);       // Rank of the sample
      nitems++;
    }
  }
  return nitems;
}

void Vecchia::_buildLHS(int nitems,
                        const std::vector<std::array<int, 4>>& neighDescr,
                        MatrixSymmetric& _matCov)
{
  SpacePoint p1;
  SpacePoint p2;

  for (int i1 = 0; i1 < nitems; i1++)
  {
    int icase1 = neighDescr[i1][0];
    int ivar1  = neighDescr[i1][1];
    int iabs1  = neighDescr[i1][3];
    if (icase1 == 1)
      _db1->getSampleAsSPInPlace(p1, iabs1);
    else
      _db2->getSampleAsSPInPlace(p1, iabs1);

    for (int i2 = 0; i2 <= i1; i2++)
    {
      int icase2 = neighDescr[i2][0];
      int ivar2  = neighDescr[i2][1];
      int iabs2  = neighDescr[i2][3];
      if (icase2 == 1)
        _db1->getSampleAsSPInPlace(p2, iabs2);
      else
        _db2->getSampleAsSPInPlace(p2, iabs2);

      double value = _model->evalCov(p1, p2, ivar1, ivar2);
      _matCov.setValue(i1, i2, value);
    }
  }
}

void Vecchia::_buildRHS(int icase2,
                        int iabs2,
                        int ivar2,
                        int nitems,
                        const std::vector<std::array<int, 4>>& neighDescr,
                        MatrixDense& _vectCov)
{
  SpacePoint p1;
  SpacePoint p2;

  if (icase2 == 1)
    _db1->getSampleAsSPInPlace(p2, iabs2);
  else
    _db2->getSampleAsSPInPlace(p2, iabs2);

  for (int i1 = 0; i1 < nitems; i1++)
  {
    int icase1 = neighDescr[i1][0];
    int ivar1  = neighDescr[i1][1];
    int iabs1  = neighDescr[i1][3];
    if (icase1 == 1)
      _db1->getSampleAsSPInPlace(p1, iabs1);
    else
      _db2->getSampleAsSPInPlace(p1, iabs1);

    double value = _model->evalCov(p1, p2, ivar1, ivar2);
    _vectCov.setValue(i1, 0, value);
  }
}

/**
 * @brief Store in the member _Y the list of data information
 * This information corresponds to the list of data defined values over the active samples
 * for all the variables
 */
void Vecchia::_loadDataFlattened()
{
  int icase = _getCase();
  int ecr   = 0;
  if (icase == 1)
  {
    _Y.resize(_Ntot1);
    int nvar = _cumulRanks1.size();
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int iech = 0, nech = _varRanks1[ivar].size(); iech < nech; iech++)
        _Y[ecr++] = _db1->getLocVariable(ELoc::Z, _varRanks1[ivar][iech], ivar);
  }
  else
  {
    _Y.resize(_Ntot2);
    int nvar = _cumulRanks2.size();
    for (int ivar = 0; ivar < nvar; ivar++)
      for (int iech = 0, nech = _varRanks2[ivar].size(); iech < nech; iech++)
        _Y[ecr++] = _db2->getLocVariable(ELoc::Z, _varRanks2[ivar][iech], ivar);
  }
}

/**
 * @brief Construct the Vecchia approximation starting from 'Ranks'
 *
 * @param Ranks MatrixT<int> which the ranks of the sample indices for each target
 * @param verbose Verbose flag
 * @return int Error returned code
 *
 * @note The dimension of 'Rank' is:
 * - ncols = Dimension of the Neighborhood + 1
 * - nrows = Number of samples (dbin and dbout [optional])
 *
 * @note For each row, the first element of 'Ranks' is the sample number
 * - if smaller than N_dbin, it refers to the sample absolute rank in 'dbin'
 * - if larger, its value (after subtracting N_dbin) gives the sample absolute
 *   rank in 'dbout'
 */
int Vecchia::computeLower(const MatrixT<int>& Ranks, bool verbose)
{
  bool debug = false;

  // Preliminary check
  int ndim;
  if (!haveSameNDim(_db1, _db2, &ndim)) return 1;
  if (ndim != (int)_model->getNDim())
  {
    messerr("Db(%d) and Model(%d) should share the same Space Dimension",
            ndim, _model->getNDim());
    return 1;
  }

  int nvar;
  if (!haveCompatibleNVar(_db1, _db2, &nvar)) return 1;
  if (nvar != (int)_model->getNVar())
  {
    messerr("Db(%d) and Model(%d) should share the same Number of Variables",
            nvar, _model->getNVar());
    return 1;
  }

  int nsample  = (int)Ranks.getNRows();
  int nb_neigh = (int)Ranks.getNCols() - 1;

  _Ntot1 = 0;
  if (_db1 != nullptr)
    _Ntot1 = _db1->getListOfSampleIndicesInPlace(nvar, _cumulRanks1, _varRanks1, true);
  _Ntot2 = 0;
  if (_db2 != nullptr)
    _Ntot2 = _db2->getListOfSampleIndicesInPlace(nvar, _cumulRanks2, _varRanks2, true);
  _Ndb1 = (_db1 != nullptr) ? _db1->getNSample() : 0;

  // Resizing
  int ntot = _Ntot1 + _Ntot2;
  _DFull.resize(ntot);
  _LFull = MatrixSparse(ntot, ntot, nb_neigh + 1);

  // Define the vector of flattened multivariate data information
  _loadDataFlattened();

  // Loop on the samples
  int nmax = nb_neigh * nvar + nvar; // Multivariate neighborhood + Collocation
  std::vector<std::array<int, 4>> neighDescr(nmax);
  for (int ivar = 0; ivar < nvar; ivar++)
  {
    for (int isample = 0; isample < nsample; isample++)
    {
      int target = Ranks(isample, 0);
      int icase1 = _getSampleCase(target);
      int iabs1  = _getAddressAbsolute(target);
      int irel1  = _getAddressInMatrix(target, ivar);

      // Build the list of neighboring information
      int nitems = _buildNeighborhood(Ranks, isample, ivar, nb_neigh, neighDescr);

      // Optional printout
      if (debug)
      {
        message("Row=%d Case=%d Variable=%d Sample=%d\n",
                irel1, icase1, ivar, iabs1);
        for (int item = 0; item < nitems; item++)
          message("- Column=%d Case=%d Variable=%d Sample=%d\n",
                  neighDescr[item][2], neighDescr[item][0],
                  neighDescr[item][1], neighDescr[item][3]);
      }

      // Fill the full matrix
      _LFull.setValue(irel1, irel1, 1.);
      double varK = _model->eval0(ivar, ivar);
      if (nitems <= 0)
      {
        // Case with no previous information available

        _DFull[irel1] = 1. / varK;
      }
      else
      {
        // Constitute the local Kriging matrix
        _vectCov.resize(nitems, 1);
        _matCov->resize(nitems, nitems);
        _buildLHS(nitems, neighDescr, *_matCov);
        _buildRHS(icase1, iabs1, ivar, nitems, neighDescr, _vectCov);

        // Solve the local system
        _chol->setMatrix(*_matCov);
        constvect vect = _vectCov.getViewOnColumn(0);
        _work.resize(vect.size());
        _chol->solve(vect, _work);

        // Patch the global matrix
        for (int i = 0; i < nitems; i++)
        {
          int irel2 = neighDescr[i][2];
          _LFull.setValue(irel2, irel1, -_work[i]);
        }
        _DFull[irel1] = 1. / (varK - VH::innerProduct(_work, vect));
      }
    }
  }

  // Final calculations
  _LFull.transposeInPlace();

  // Optional printout
  if (verbose)
  {
    mestitle(1, "Lower Vecchia Matrix");
    _LFull.display();
    VH::dump("Diagonal of Vecchia Matrix", _DFull);
  }
  return 0;
}

// Calculate LdY = Ldat %*% Y
VectorDouble Vecchia::calculateLdY(const VectorDouble& Y) const
{
  int nd = getND();
  int nt = getNT();
  VectorDouble LdY(nd);

  for (int id = 0; id < nd; id++)
  {
    double value = 0.;
    for (int jd = 0; jd < nd; jd++)
      value += getLFull(id + nt, jd + nt) * Y[jd];
    LdY[id] = value;
  }
  return LdY;
}

// Calculate FtLdY = Ft %*% Ldat %*% Y
VectorDouble Vecchia::calculateFtLdY(const VectorDouble& LdY) const
{
  int nd = getND();
  int nt = getNT();
  VectorDouble FtLdY(nt);
  for (int it = 0; it < nt; it++)
  {
    double value = 0.;
    for (int id = 0; id < nd; id++)
      value += getLFull(id + nt, it) * LdY[id];
    FtLdY[it] = value;
  }
  return FtLdY;
}

MatrixSparse* Vecchia::calculateW(const VectorDouble& D_dd) const
{
  int nd = getND();
  int nt = getNT();

  // Extract sub-part of 'Diagonal' vector
  VectorDouble D_tt(nt);
  VH::extractInPlace(getDFull(), D_tt, 0);

  // Extracting information from 'LFull'
  VectorInt indT = VectorInt(nt + nd, -1);
  for (int it = 0; it < nt; it++) indT[it] = it;
  VectorInt indD = VectorInt(nt + nd, -1);
  for (int id = 0; id < nd; id++) indD[id + nt] = id;
  MatrixSparse* Ltt = getLFull().extractSubmatrixByRanks(indT, indT);
  Ltt->forceDimension(nt, nt);
  MatrixSparse* Ldt = getLFull().extractSubmatrixByRanks(indD, indT);
  Ldt->forceDimension(nd, nt);

  /*! Product 't(A)' %*% 'M' %*% 'A' or 'A' %*% 'M' %*% 't(A)' */
  MatrixSparse* mat1 = prodNormMat(Ltt, D_tt, true);
  MatrixSparse* mat2 = prodNormMat(Ldt, D_dd, true);
  mat1->forceDimension(nt, nt);
  mat2->forceDimension(nt, nt);

  // Cleaning
  indT.clear();
  indD.clear();
  delete Ltt;
  delete Ldt;

  MatrixSparse* W = MatrixSparse::addMatMat(mat1, mat2);
  delete mat1;
  delete mat2;
  return W;
}

int krigingVecchia(Db* dbin,
                   Db* dbout,
                   ModelGeneric* model,
                   int nb_neigh,
                   bool verbose,
                   const NamingConvention& namconv)
{
  Vecchia V = Vecchia(model, nb_neigh, dbout, dbin);

  MatrixT<int> Ranks = findNN(dbout, dbin, nb_neigh + 1, false, verbose);
  if (V.computeLower(Ranks, verbose)) return 1;

  // Extract sub-part of 'Diagonal' vector
  VectorDouble DFull = V.getDFull();
  int nd             = V.getND();
  int nt             = V.getNT();
  int nvar           = model->getNVar();
  VectorDouble D_dd(nd);
  VH::extractInPlace(DFull, D_dd, nt);

  // Calculate LdY
  const VectorDouble& Y = V.getY();
  VectorDouble LdY      = V.calculateLdY(Y);
  VH::multiplyInPlace(LdY, D_dd);

  // Calculate FtLdY
  VectorDouble FtLdY = V.calculateFtLdY(LdY);

  // Calculating 'W'
  MatrixSparse* W = V.calculateW(D_dd);

  // Compute the Cholesky decomposition of 'W'
  CholeskySparse cholW(*W);
  if (!cholW.isReady())
  {
    messerr("Cholesky decomposition of Covariance matrix failed");
    return 1;
  }

  // Perform the estimation
  VectorDouble result = cholW.solveX(FtLdY);
  for (int i = 0; i < nt; i++) result[i] = -result[i];

  // Saving the results
  int iptr = dbout->addColumns(result, String(), ELoc::UNKNOWN, 0, true);
  namconv.setNamesAndLocators(dbout, iptr, "estim", nvar);

  return 0;
}

void Vecchia::productVecchia(constvect Y, vect res) const
{
  _LdY.resize(_LFull.getNRows());
  _LFull.prodMatVecInPlace(Y, _LdY, false);
  VH::multiplyInPlace(_LdY, _DFull);
  _LFull.prodMatVecInPlace(_LdY, res, true);
}

void Vecchia::productMatVecchia(const MatrixDense& X, MatrixDense& resmat) const
{
  int nrows = X.getNRows();
  int ncols = X.getNCols();
  resmat.resize(nrows, ncols);

  // Loop on the columns
  for (int icol = 0; icol < ncols; icol++)
  {
    constvect colin = X.getViewOnColumn(icol);
    vect colout     = resmat.getViewOnColumnModify(icol);
    productVecchia(colin, colout);
  }
}

/**
 * Compute the log-likelihood (based on Vecchia approximation for covMat)
 *
 * @param db  Db structure where variable are loaded from
 * @param model ModelGeneric structure used for the calculation
 * @param nb_neigh Number of neighbors to consider in the Vecchia approximation
 * @param verbose Verbose flag
 *
 * @remarks The calculation considers all the active samples.
 * @remarks It can work in multivariate case with or without drift conditions (linked or not)
 * @remarks The algorithm is stopped (with a message) in the heterotopic case
 */
double logLikelihoodVecchia(const Db* db,
                            ModelGeneric* model,
                            int nb_neigh,
                            bool verbose)
{
  Vecchia* vec  = Vecchia::createForOptim(model, db, nb_neigh);
  double result = vec->computeCost(verbose);
  delete vec;
  return result;
}

Vecchia* Vecchia::createForOptim(ModelGeneric* model,
                                 const Db* db,
                                 int nb_neigh,
                                 bool reml)
{

  auto* vec            = new Vecchia(model, nb_neigh, db, nullptr, reml);
  MatrixSymmetric vars = dbVarianceMatrix(db);
  double hmax          = db->getExtensionDiagonal();
  vec->setEnvironment(vars, hmax);
  vec->init();
  return vec;
}

void Vecchia::_computeCm1X()
{
  productMatVecchia(_X, _Cm1X);
}

void Vecchia::_computeCm1Y()
{
  _Cm1Y.resize(_Y.size());
  productVecchia(_Y, _Cm1Y);
}

double Vecchia::_computeLogDet() const
{
  return -VH::cumulLog(getDFull());
}

void Vecchia::_updateModel(bool verbose)
{
  computeLower(_Ranks, verbose);
}
} // namespace gstlrn