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
#include "geoslib_define.h"

#include "Enum/EAnam.hpp"

#include "Anamorphosis/AnamDiscrete.hpp"
#include "Matrix/MatrixSquare.hpp"
#include "Stats/PCA.hpp"
#include "Stats/Selectivity.hpp"

class GSTLEARN_EXPORT AnamDiscreteDD: public AnamDiscrete
{
public:
  AnamDiscreteDD(double mu = 1., double scoef = 0.);
  AnamDiscreteDD(const AnamDiscreteDD &m);
  AnamDiscreteDD& operator= (const AnamDiscreteDD &m);
  virtual ~AnamDiscreteDD();

  /// ICloneable Interface
  IMPLEMENT_CLONING(AnamDiscreteDD)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASerializable Interface
  static AnamDiscreteDD* createFromNF(const String& neutralFilename, bool verbose = true);

  /// AAnam Interface
  const EAnam&  getType() const override { return EAnam::fromKey("DISCRETE_DD"); }
  bool hasFactor() const override { return true; }
  int getNFactor() const override { return 0; }
  VectorDouble z2factor(double z, const VectorInt& ifacs) const override;
  double computeVariance(double sval) const override;
  int  updatePointToBlock(double r_coef) override;
  bool allowChangeSupport() const override { return true; }
  bool isChangeSupportDefined() const override { return (_sCoef > 0.); }
  int fitFromArray(const VectorDouble &tab,
                   const VectorDouble &wt = VectorDouble()) override;

  /// AnamDiscrete Interface
  void calculateMeanAndVariance() override;

  VectorDouble factors_exp(bool verbose = false);
  VectorDouble factors_maf(bool verbose = false);
  VectorDouble factors_mod();
  MatrixSquare chi2I(const VectorDouble& chi, int mode);

  static AnamDiscreteDD* create(double mu = 1., double scoef = 0.);
  void reset(int ncut,
             double scoef,
             double mu,
             const VectorDouble &zcut,
             const MatrixSquare &pcaz2f,
             const MatrixSquare &pcaf2z,
             const VectorDouble &stats);

  PCA& getMAF() { return _maf; }
  double getMu() const { return _mu; }
  double getSCoef() const { return _sCoef; }
  const MatrixSquare& getI2Chi() const { return _i2Chi; }
  MatrixSquare getPcaZ2Fs() const { return _maf.getZ2Fs(); }
  MatrixSquare getPcaF2Zs() const { return _maf.getF2Zs(); }

  void setMu(double mu) { _mu = mu; }
  void setRCoef(double rcoef) { _sCoef = rcoef; }
  void setPcaZ2F(const MatrixSquare& pcaz2f) { _maf.setZ2Fs(pcaz2f); }
  void setPcaF2Z(const MatrixSquare& pcaf2z) { _maf.setF2Zs(pcaf2z); }
  void setI2Chi(const MatrixSquare& i2Chi) { _i2Chi = i2Chi; }

  int factor2Selectivity(Db *db,
                         Selectivity* selectivity,
                         const VectorInt& cols_est,
                         const VectorInt& cols_std,
                         int iptr0);

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "AnamDiscreteDD"; }

private:
  int _stats(int nech, const VectorDouble& tab);
  VectorDouble _generator(const VectorDouble& vecc,
                          const VectorDouble& veca,
                          const VectorDouble& vecb,
                          VectorDouble& eigvec,
                          VectorDouble& eigval);
  void _lambdaToMul();
  void _blockAnamorphosis(const VectorDouble& chi);
  void _globalSelectivity(Selectivity* selectivity);

private:
  double _mu;
  double _sCoef;
  PCA    _maf;
  MatrixSquare _i2Chi; // Dimension: nclass * nfacies

  friend class Selectivity;
};
