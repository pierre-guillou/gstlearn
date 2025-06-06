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

#include "Enum/EDirGen.hpp"
#include "Enum/EGaussInv.hpp"

#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/NamingConvention.hpp"

class Db;
class MatrixSymmetric;
class MatrixDense;
class AMatrix;

class GSTLEARN_EXPORT PPMT : public AStringable, public ICloneable
{
public:
  PPMT(int ndir = 50,
       bool flagPreprocessing = false,
       const EDirGen& methodDir = EDirGen::fromKey("VDC"),
       const EGaussInv& methodTrans = EGaussInv::fromKey("EMP"),
       int nbpoly = 30,
       double alpha = 2.);
  PPMT(const PPMT &m);
  PPMT& operator= (const PPMT &m);
  virtual ~PPMT();

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ICloneable Interface
  IMPLEMENT_CLONING(PPMT)

  static PPMT* create(int ndir = 50,
                      bool flagPreprocessing = false,
                      const EDirGen& methodDir = EDirGen::fromKey("VDC"),
                      const EGaussInv& methodTrans = EGaussInv::fromKey("EMP"),
                      int nbpoly = 30,
                      double alpha = 2.);

  int getNiter()       const { return _niter;  }
  double getAlpha()    const { return _alpha;  }
  int getNdir()        const { return _ndir;   }
  int getNdim()        const { return _ndim;   }
  int getNbpoly()      const { return _nbpoly; }
  bool isFitted()      const { return _isFitted; }
  const EDirGen& getMethodDir() const { return _methodDir; }
  const EGaussInv& getMethodTrans() const { return _methodTrans; }

  VectorDouble getSerieAngle() const { return _serieAngle; }
  VectorDouble getSerieScore(bool flagLog = false) const;

  int fit(Db *db,
          const VectorString &names,
          bool flagStoreInDb = false,
          int niter = 100,
          bool verbose = false,
          const NamingConvention &namconv = NamingConvention("Y"));
  int fitFromMatrix(AMatrix *Y, int niter, bool verbose = false);
  int rawToGaussian(Db *db,
                    const VectorString &names,
                    int niter = 0,
                    const NamingConvention &namconv = NamingConvention("Y"));
  int gaussianToRaw(Db *db,
                    const VectorString &names,
                    int niter = 0,
                    const NamingConvention &namconv = NamingConvention("Z"));

private:
  void _generateAllDirections();
  void _fitInitHermite(AMatrix* Y);
  void _initGaussianizeForward(AMatrix* Y);
  void _initGaussianizeBackward(AMatrix* Y);
  void _iterationFit(AMatrix *Y, const VectorDouble &N0);
  void _iterationForward(AMatrix *Y, const VectorDouble &N0, int iter = 0);
  void _iterationBackward(AMatrix *Y, const VectorDouble &N0, int iter = 0);
  static double _gaussianizeForward(double Yi,
                                    int rank,
                                    const AnamHermite* anam,
                                    const VectorDouble& N0);
  static double _gaussianizeBackward(double Yi, const AnamHermite *anam);
  void _projectOnDirection(const AMatrix *Y, int id, VectorDouble &Y0);
  double _getGaussianDistance(const VectorDouble &Yi,
                              const VectorInt &Ri,
                              const VectorDouble &N0) const;
  void _shiftForward(AMatrix *Y,
                     int id,
                     const AnamHermite *anam,
                     const VectorDouble &Y0,
                     const VectorInt &R0,
                     const VectorDouble &N0) const;
  void _shiftBackward(AMatrix *Y,
                      int id,
                      const AnamHermite *anam,
                      const VectorDouble &Y0) const;

private:
  int _niter;
  int _ndir;
  int _nbpoly;
  double _alpha;
  EDirGen _methodDir;
  EGaussInv _methodTrans;
  bool _flagPreprocessing;

  bool _isFitted;

  mutable int _ndim;
  mutable VectorDouble _serieAngle;
  mutable VectorDouble _serieScore;
  mutable MatrixDense*   _dirmat;
  mutable std::vector<AnamHermite*> _anams;
  mutable std::vector<AnamHermite*> _initAnams;
  mutable MatrixDense* _initSphering;
};

