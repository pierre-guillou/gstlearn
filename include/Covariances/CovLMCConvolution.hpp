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

#include "Covariances/CovContext.hpp"
#include "gstlearn_export.hpp"

#include "Enum/EConvDir.hpp"
#include "Enum/EConvType.hpp"

#include "Covariances/CovAnisoList.hpp"
#include "Matrix/MatrixDense.hpp"

#include "Space/ASpace.hpp"

class SpacePoint;
class CovAniso;
class Model;

typedef struct
{
  String convName;
  double convScale;
  double (*convFunc)(double);
} Def_Convolution;

/* Prototyping the internal covariance functions */
GSTLEARN_EXPORT double _conv_uniform(double v);
GSTLEARN_EXPORT double _conv_exponential(double v);
GSTLEARN_EXPORT double _conv_gaussian(double v);
GSTLEARN_EXPORT double _conv_sincard(double v);

GSTLEARN_EXPORT Def_Convolution& D_CONV(int rank);

class GSTLEARN_EXPORT CovLMCConvolution : public CovAnisoList
{
public:
  CovLMCConvolution(const EConvType& conv_type,
                    const EConvDir& conv_dir,
                    double conv_range,
                    int conv_ndisc,
                    const CovContext& ctxt = CovContext());
  CovLMCConvolution(const CovLMCConvolution &r);
  CovLMCConvolution& operator= (const CovLMCConvolution &r);
  virtual ~CovLMCConvolution();

  /// ICloneable interface
  IMPLEMENT_CLONING(CovLMCConvolution)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ACov Interface
  virtual double eval0(int ivar = 0,
                       int jvar = 0,
                       const CovCalcMode* mode = nullptr) const override;
  
  
  int init(const EConvType& conv_type, const EConvDir& conv_idir, double conv_range, int conv_ndisc);

  double getConvRange() const { return _convRange; }
  const VectorDouble& getConvWeight() const { return _convWeight; }
  const MatrixDense& getConvIncr() const { return _convIncr; }
  VectorDouble getConvIncr(int rank) const { return _convIncr.getColumn(rank); }
  int getNConv() const { return _convNumber; }

protected:
  void _optimizationSetTarget(SpacePoint& pt) const override
  {
    ACov::_optimizationSetTarget(pt); 
  }

private:
virtual double _eval(const SpacePoint& p1,
                     const SpacePoint& p2,
                     int ivar = 0,
                     int jvar = 0,
                     const CovCalcMode* mode = nullptr) const override;
private:
  EConvType _convType; /* Convolution type */
  EConvDir  _convDir;  /* Convolution direction: 0:X, 1:Y, 2:Z, 3:XY, 4:XYZ */
  int _convDiscNumber; /* Number of discretization per direction */
  double _convRange; /* Convolution Range */
  int _convNumber;
  MatrixDense _convIncr; /* Discretization lags */
  VectorDouble      _convWeight; /* Weights for convolution */
};
