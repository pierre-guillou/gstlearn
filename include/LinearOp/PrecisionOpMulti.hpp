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

#include "LinearOp/ALinearOp.hpp"
#include "LinearOp/ASimulable.hpp"
#include "gstlearn_export.hpp"

#include "Basic/VectorT.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "Model/Model.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"
#include <vector>


#define IND(i,j,nvar) j * nvar + i - (j * (j + 1))/2

class Model;
class ASimulable;

/**
 * Class to store objects for SPDE
 */
class GSTLEARN_EXPORT PrecisionOpMulti : public AStringable, public ASimulable
{
public:
  PrecisionOpMulti(Model* model               = nullptr,
                   const VectorMeshes& meshes = VectorMeshes(),
                   bool buildOp               = true);
  PrecisionOpMulti(const PrecisionOpMulti& m)            = delete;
  PrecisionOpMulti& operator=(const PrecisionOpMulti& m) = delete;
  virtual ~PrecisionOpMulti();
  MatrixSquareSymmetric getInvCholSill(int icov) const
  {
    return _invCholSills[icov];
  }
  int getSize() const override;

protected:
  void buildQop();
#ifndef SWIG

protected:
  int _addToDest(const constvect vecin, vect vecout) const override;
  int _addSimulateToDest(const constvect vecin, vect vecout) const override;
#endif

public:
  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  VectorDouble evalSimulate(const VectorDouble& vec);

protected:
  int size(int imesh) const;
  int _getNCov() const;
  int _getCovInd(int i) const { return _covList[i]; }
  int _getNVar() const;
  int _getNMesh() const;

private:
  virtual int _addToDestImpl(const constvect vecin, vect vecout) const;
  bool _checkReady() const;
  virtual void _buildQop();
  bool _isValidModel(Model* model);
  bool _isValidMeshes(const std::vector<const AMesh*>& meshes);
  bool _isNoStat(int istruct) const { return _isNoStatForVariance[istruct]; }
  bool _matchModelAndMeshes() const;

  int _buildGlobalMatricesStationary(int icov);
  int _buildLocalMatricesNoStat(int icov);
  int _buildMatrices();
  void _popsClear();
  void _computeSize();

protected:
  std::vector<PrecisionOp*> _pops;
  VectorBool _isNoStatForVariance;
  std::vector<VectorVectorDouble> _invCholSillsNoStat;
  std::vector<VectorVectorDouble> _cholSillsNoStat;
  std::vector<MatrixSquareSymmetric> _invCholSills; // Inverse of Cholesky of the Sills
  std::vector<MatrixSquareSymmetric> _cholSills; // Cholesky of the Sills
  Model* _model; // Not to be deleted. TODO : make it const
  std::vector<const AMesh*> _meshes; // Not to be deleted
  int _size;

private:
  bool _isValid;
  VectorInt _covList;
  VectorInt _nmeshList;
  bool _allStat;
  bool _ready;

private:
  mutable VectorVectorDouble _works;
  mutable VectorDouble _workTot;
};
