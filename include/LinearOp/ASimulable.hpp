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

#include "Basic/VectorNumT.hpp"
#include "Matrix/VectorEigen.hpp"
#include "LinearOp/ALinearOp.hpp"

#ifndef SWIG
  #include <Eigen/Core>
  #include <Eigen/Dense>
#endif

class GSTLEARN_EXPORT ASimulable : public ALinearOp
{
public:
  ASimulable() { }
  virtual ~ASimulable() {}
  
  int evalSimulate(const VectorDouble& whitenoise, VectorDouble& outv) const;
  VectorDouble evalSimulate(const VectorDouble& whitenoise) const;
  int evalSimulate(const VectorEigen& whitenoise, VectorEigen& outv) const;
  int addSimulateToDest(const VectorDouble& whitenoise,
                              VectorDouble& outv) const;
  int addSimulateToDest(const VectorEigen& whitenoise, VectorEigen& outv) const;
#ifndef SWIG

public:
  int evalSimulate(const Eigen::VectorXd& whitenoise,
                         Eigen::VectorXd& outv) const;
  int addSimulateToDest(const Eigen::VectorXd& whitenoise,
                              Eigen::VectorXd& outv) const;

protected:
  virtual int _addSimulateToDest(const Eigen::VectorXd& whitenoise,
                                       Eigen::VectorXd& outv) const = 0;
#endif
};