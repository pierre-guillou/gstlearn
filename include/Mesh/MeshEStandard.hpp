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

#include "Basic/VectorNumT.hpp"
#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixInt.hpp"

class MeshETurbo;

/**
 * Standard Meshing defined in the Euclidean space
 */
class GSTLEARN_EXPORT MeshEStandard: public AMesh
{
public:
  MeshEStandard();
  MeshEStandard(const MeshEStandard &m);
  MeshEStandard& operator=(const MeshEStandard &m);
  virtual ~MeshEStandard();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for AMesh
  int     getNApices() const override;
  int     getNMeshes() const override;
  int     getApex(int imesh, int rank) const override;
  double  getCoor(int imesh, int rank, int idim) const override;
  double  getApexCoor(int i, int idim) const override;
  double  getMeshSize(int imesh) const override;
  static MeshEStandard* createFromNF(const String& neutralFilename,
                                     bool verbose = true);
  static MeshEStandard* createFromExternal(const MatrixDense& apices,
                                           const MatrixInt& meshes,
                                           bool verbose = false);

  VectorInt    getMeshList() const { return _meshes.getValues(); }
  VectorDouble getPointList(bool byCol = true) const;
  int reset(const MatrixDense& apices,
            const MatrixInt& meshes,
            bool verbose = false);
  int reset(int ndim,
            int napexpermesh,
            const VectorDouble &apices,
            const VectorInt &meshes,
            bool byCol = true,
            bool verbose = false);
  int resetFromTurbo(const MeshETurbo &turbo, bool verbose = false);

  const MatrixDense& getApices() const { return _apices; }
  const MatrixInt& getMeshes() const { return _meshes; }

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os,bool verbose = false) const override;
  String _getNFName() const override { return "MeshEStandard"; }
  void _defineBoundingBox(void);

private:
  bool _coorInMesh(const VectorDouble& coor,
                   int imesh,
                   double meshsize,
                   VectorDouble& weights) const;
  void _deallocate();
  int  _recopy(const MeshEStandard &m);
  void _checkConsistency() const;
  void _setApex(int imesh, int rank, int value);
  void _validate();

private:
  MatrixDense _apices; // Dimension: NRow=napices; Ncol=Ndim
  MatrixInt         _meshes; // Dimension: Nrow=Nmesh; Ncol=NApexPerMesh
};
