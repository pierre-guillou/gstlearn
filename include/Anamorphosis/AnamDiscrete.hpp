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

#include "Anamorphosis/AAnam.hpp"
#include "Matrix/MatrixDense.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"

namespace gstlrn
{
class Selectivity;

class AAnam; // Forward declaration

class GSTLEARN_EXPORT AnamDiscrete: public AAnam
{
public:
  AnamDiscrete();
  AnamDiscrete(const AnamDiscrete& m);
  AnamDiscrete& operator=(const AnamDiscrete& m);
  virtual ~AnamDiscrete();

  /// AStringable Interface
  String toString(const AStringFormat* strfmt = nullptr) const override;

  /// AAnam interface
  bool hasGaussian() const override { return false; }
  int getNClass() const override { return _nCut + 1; }

  /// Interface for AnamDiscrete
  virtual void calculateMeanAndVariance();
  double getVariance() const override { return _variance; }

  int getNCut() const { return _nCut; }
  int getNElem() const { return _nElem; }
  const VectorDouble& getZCut() const { return _zCut; }
  double getZCut(int i) const { return _zCut[i]; }
  double getMean() const { return _mean; }

  void setMean(double mean) { _mean = mean; }
  void setVariance(double variance) { _variance = variance; }
  void setNCut(int ncut);
  void setZCut(const VectorDouble& zcut);
  void setNElem(int nelem);
  void setStats(const VectorDouble& stats);

  // Function for using Stats in DD anamorphosis
  double getDDStatProp(int iclass) const;
  double getDDStatZmoy(int iclass) const;
  double getDDStatCnorm(int iclass) const;
  double getDDStatLambda(int iclass) const;
  double getDDStatU(int iclass) const;
  double getDDStatMul(int iclass) const;
  void setDDStatProp(int iclass, double value);
  void setDDStatZmoy(int iclass, double value);
  void setDDStatCnorm(int iclass, double value);
  void setDDStatLambda(int iclass, double value);
  void setDDStatU(int iclass, double value);
  void setDDStatMul(int iclass, double value);

  // Function for using Stats in IR anamorphosis
  double getIRStatT(int iclass) const;
  double getIRStatQ(int iclass) const;
  double getIRStatZ(int iclass) const;
  double getIRStatB(int iclass) const;
  double getIRStatR(int iclass) const;
  double getIRStatRV(int iclass) const;
  void setIRStatT(int iclass, double value);
  void setIRStatQ(int iclass, double value);
  void setIRStatZ(int iclass, double value);
  void setIRStatB(int iclass, double value);
  void setIRStatR(int iclass, double value);
  void setIRStatRV(int iclass, double value);

  const MatrixDense& getStats() const { return _stats; }

protected:
  bool _deserializeAscii(std::istream& is, bool verbose = false) override;
  bool _serializeAscii(std::ostream& os, bool verbose = false) const override;
#ifdef HDF5
  bool _deserializeH5(H5::Group& grp, bool verbose = false) override;
  bool _serializeH5(H5::Group& grp, bool verbose = false) const override;
#endif
  String _getNFName() const override { return "AnamDiscrete"; }

  bool _isClassValid(int iclass) const;
  void _resize();

private:
  int _nCut;
  int _nElem;
  double _mean;
  double _variance;
  VectorDouble _zCut;
  MatrixDense _stats;
};
} // namespace gstlrn