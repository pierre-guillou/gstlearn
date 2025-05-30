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

#include "Space/ASpace.hpp"
#include "Basic/VectorNumT.hpp"
#include <memory>

class SpacePoint;

class GSTLEARN_EXPORT SpaceSN: public ASpace
{
private:
  SpaceSN(unsigned int ndim, double radius);
  SpaceSN(const SpaceSN &r);
  SpaceSN& operator=(const SpaceSN &r);
 
public:
  virtual ~SpaceSN();

  /// ICloneable interface
  IMPLEMENT_CLONING(SpaceSN)

  /// Return the concrete space type
  ESpaceType getType() const override { return ESpaceType::SN; };

  static ASpaceSharedPtr create(int ndim, double radius);
  /// Return the sphere radius
  double getRadius() const { return _radius; }

  /// Dump a space in a string
  virtual String toString(const AStringFormat* strfmt,
                          int idx = -1) const override;

  /// Return true if the given space is equal to me
  virtual bool isEqual(const ASpace *space) const override;

protected:

  /// Move the given space point by the given vector
  void _move(SpacePoint &p1, const VectorDouble &vec) const override;

  /// Return the distance between two space points
  double _getDistance(const SpacePoint& p1,
                      const SpacePoint& p2,
                      int ispace = -1) const override;

  /// Return the distance between two space points with the given tensor
  double _getDistance(const SpacePoint& p1,
                      const SpacePoint& p2,
                      const Tensor& tensor,
                      int ispace = -1) const override;

  /// Return the distance in frequential domain between two space points with the given tensor
  double _getFrequentialDistance(const SpacePoint& p1,
                                 const SpacePoint& p2,
                                 const Tensor& tensor,
                                 int ispace = -1) const override;

  /// Return the increment vector between two space points for the current space context
  VectorDouble _getIncrement(const SpacePoint& p1,
                             const SpacePoint& p2,
                             int ispace = -1) const override;
  
  /// Return the increment vector between two space points in a given vector
  void _getIncrementInPlace(const SpacePoint& p1,
                            const SpacePoint& p2,
                            VectorDouble& ptemp,
                            int ispace = -1) const override;

private:
  /// Sphere radius
  double _radius;
};

