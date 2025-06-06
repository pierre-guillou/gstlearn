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

#include "Enum/EPowerPT.hpp"
#include "Basic/ICloneable.hpp"
#include "Mesh/AMesh.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"

#include <memory>

#ifndef SWIG

#include <Eigen/Core>
#include <Eigen/Dense>
#endif

class CovAniso;
class EConsElem;

/**
 * \brief Shift Operator for performing the basic tasks of SPDE
 */

#ifndef SWIG
#  include "LinearOp/ALinearOpEigenCG.hpp"
DECLARE_EIGEN_TRAITS(AShiftOp)
#else
#  include "LinearOp/ALinearOp.hpp"
#endif

class GSTLEARN_EXPORT AShiftOp: public ICloneable,
#ifndef SWIG
  public ALinearOpEigenCG<AShiftOp>
#else
  public ALinearOp
#endif
{
public:
  AShiftOp(CovAniso* cova = nullptr, int napices = 0);
  AShiftOp(const AShiftOp& shift);
  AShiftOp& operator=(const AShiftOp& shift);
  virtual void prodLambda(const VectorDouble& x,
                          VectorDouble& y,
                          const EPowerPT& power) const;
  virtual ~AShiftOp();
  virtual double getMaxEigenValue() const = 0;

  virtual void normalizeLambdaBySills(const AMesh*) = 0;
  const VectorDouble& getLambdas() const { return _Lambda; }
  virtual double getLambda(int iapex) const { return _Lambda[iapex]; }

  static std::shared_ptr<CovAniso> cloneAndCast(const CovAniso* cova);
  static std::shared_ptr<CovAniso> cloneAndCast(const std::shared_ptr<CovAniso> &cova);
  int getSize() const override { return _napices; }

#ifndef SWIG
    virtual void addProdLambda(const constvect x, vect y, const EPowerPT& power) const;
    void prodLambda(const constvect x, vect y, const EPowerPT& power) const;
    void prodLambda(const VectorDouble& x, vect y, const EPowerPT& power) const;
    void prodLambda(const constvect x, VectorDouble& y, const EPowerPT& power) const;
#endif
#ifndef SWIG
    int _addToDest(const constvect inv, vect outv) const override = 0;
#endif

protected:
    std::shared_ptr<CovAniso>& _getCovAniso();
    void _setCovAniso(const CovAniso* cova);
    bool _isNoStat();
    bool _isGlobalHH();

protected:
    VectorDouble _Lambda;
    int _napices;
    // Following list of members are there to ease the manipulation and reduce
    // argument list
    std::shared_ptr<CovAniso> _cova;
};

#ifndef SWIG
  DECLARE_EIGEN_PRODUCT(AShiftOp)
#endif
