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

#include "Basic/NamingConvention.hpp"
#include "Basic/VectorNumT.hpp"
#include "Estimation/ALikelihood.hpp"
#include "LinearOp/CholeskyDense.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixSparse.hpp"
#include "Matrix/MatrixSymmetric.hpp"
#include "Matrix/MatrixT.hpp"
#include "Mesh/AMesh.hpp"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class Db;
class ModelGeneric;

class GSTLEARN_EXPORT Vecchia: public ALikelihood
{
public:
  Vecchia(ModelGeneric* model,
          int nb_neigh,
          const Db* db1,
          const Db* db2 = nullptr,
          bool reml     = false);
  Vecchia(const Vecchia& r);
  Vecchia& operator=(const Vecchia& r);
  virtual ~Vecchia();

public:
  static Vecchia* createForOptim(ModelGeneric* model,
                                 const Db* db1,
                                 int nb_neigh = 30,
                                 bool reml    = false);

  int computeLower(const MatrixT<int>& Ranks, bool verbose = false);
  const MatrixSparse& getLFull() const { return _LFull; }
  const VectorDouble& getDFull() const { return _DFull; }
  const VectorDouble& getY() const { return _Y; }

  double getLFull(int i, int j) const { return _LFull.getValue(i, j); }
  int getND() const { return _Ntot2; }
  int getNT() const { return _Ntot1; }
  int getNonZeros() const { return _LFull.getNonZeros(); }

  void productMatVecchia(const MatrixDense& X, MatrixDense& resmat) const;
  void productVecchia(constvect Y, vect res) const;
  VectorDouble calculateLdY(const VectorDouble& Y) const;
  VectorDouble calculateFtLdY(const VectorDouble& LdY) const;
  MatrixSparse* calculateW(const VectorDouble& D_dd) const;

private:
  void _init(bool verbose = false) override;
  void _updateModel(bool verbose = false) override;
  void _computeCm1X() override;
  void _computeCm1Y() override;
  double _computeLogDet() const override;
  int _buildNeighborhood(const MatrixT<int>& Ranks,
                         int isample,
                         int ivar,
                         int nb_neigh,
                         std::vector<std::array<int, 4>>& neighDescr) const;
  void _buildLHS(int nitems,
                 const std::vector<std::array<int, 4>>& neighDescr,
                 MatrixSymmetric& _matCov);
  void _buildRHS(int icase2,
                 int iabs2,
                 int ivar2,
                 int nitems,
                 const std::vector<std::array<int, 4>>& neighDescr,
                 MatrixDense& _vectCov);
  void _loadDataFlattened();
  int _getAddressInMatrix(int ip, int ivar) const;
  int _getAddressAbsolute(int ip) const;
  int _getSampleCase(int ip) const;
  int _getCase() const;

private:
  // Following members are copies of pointers (not to be deleted)
  int _nbNeigh;
  const Db* _db1;
  const Db* _db2;

  MatrixT<int> _Ranks; // Matrix of ranks for the Vecchia approximation
  std::shared_ptr<MatrixSymmetric> _matCov;
  MatrixDense _vectCov;
  VectorDouble _work; // Work vector for calculations
  mutable VectorDouble _Y;
  mutable VectorDouble _LdY;
  mutable VectorDouble _DFull;
  mutable MatrixSparse _LFull;
  mutable CholeskyDense* _chol; // Cholesky decomposition of the covariance matrix
  // Local calculation results (to be deleted later)
  mutable int _Ndb1; // Number of samples in Db1 (used for shift calculations)
  mutable int _Ntot1;
  mutable int _Ntot2;
  mutable VectorInt _cumulRanks1;
  mutable VectorInt _cumulRanks2;
  mutable VectorVectorInt _varRanks1;
  mutable VectorVectorInt _varRanks2;
};

GSTLEARN_EXPORT int krigingVecchia(Db* dbin,
                                   Db* dbout,
                                   ModelGeneric* model,
                                   int nb_neigh                    = 5,
                                   bool verbose                    = false,
                                   const NamingConvention& namconv = NamingConvention("Vecchia"));
GSTLEARN_EXPORT double logLikelihoodVecchia(const Db* db,
                                            ModelGeneric* model,
                                            int nb_neigh = 5,
                                            bool verbose = false);
} // namespace gstlrn