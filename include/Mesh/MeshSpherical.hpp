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

#include "Mesh/AMesh.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Matrix/MatrixInt.hpp"

/**
 * Meshing defined in the Spherical Space
 */
class GSTLEARN_EXPORT MeshSpherical : public AMesh
{
public:
  MeshSpherical(const MatrixDense& apices = MatrixDense(),
                const MatrixInt& meshes = MatrixInt());
  MeshSpherical(const MeshSpherical &m);
  MeshSpherical& operator= (const MeshSpherical &m);
  virtual ~MeshSpherical();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface to AMesh
  int     getNApices() const override;
  int     getNMeshes() const override;
  double  getMeshSize(int imesh) const override;
  int     getApex(int imesh, int rank) const override;
  double  getCoor(int imesh, int rank, int idim) const override;
  double  getApexCoor(int i, int idim) const override;
  int     getEmbeddedNDim() const override { return 3; }
  void    getEmbeddedCoorPerMesh(int imesh, int ic, VectorDouble& coords) const override;
  void    getEmbeddedCoorPerApex(int iapex, VectorDouble& coords) const override;
  void    getBarycenterInPlace(int imesh, VectorDouble& coord) const override;

  static MeshSpherical* createFromNF(const String& neutralFilename,
                                     bool verbose = true);
  static MeshSpherical* create(const MatrixDense& apices = MatrixDense(),
                               const MatrixInt& meshes         = MatrixInt());

  int reset(int ndim,
            int napexpermesh,
            const VectorDouble &apices,
            const VectorInt &meshes,
            bool byCol,
            bool verbose = false);
  int  getVariety() const override { return 1; }

  const MatrixDense& getApices() const { return _apices; }
  const MatrixInt& getMeshes() const { return _meshes; }
  VectorVectorInt getMeshesAsVVI() const {return _meshes.getMatrix();}

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os,bool verbose = false) const override;
  String _getNFName() const override { return "MeshSpherical"; }

private:
  void _defineBoundingBox();
  VectorDouble _defineUnits() const;
  int _recopy(const MeshSpherical &m);
  static double _closestValue(double ref, double coor, double period);
  void _checkConsistency() const;
  bool _weightsInMesh(const VectorDouble& coor,
                      const VectorVectorDouble& corners,
                      double meshsize,
                      VectorDouble& weights,
                      double eps = EPSILON5) const override;
  static void _getCoordOnSphere(double longitude,
                                double latitude,
                                VectorDouble& coords);

private:
  MatrixDense _apices; // Dimension: NRow=napices; Ncol=Ndim(=2)
  MatrixInt         _meshes; // Dimension: Nrow=Nmesh; Ncol=NApexPerMesh
};
