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

#include "Space/ASpace.hpp"
#include "gstlearn_export.hpp"

#include "Enum/ESpaceType.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/VectorNumT.hpp"

class ASpace;
class SpacePoint;

/**
 * This class is the base class for all objects that need to know what is its space definition.
 * All ASpaceObject can access to the number of space dimensions and can ask to calculate
 * a distance between two ASpaceObjects.
 *
 * This class also stores a unique (static) default global space that will be used as default space
 * when creating a ASpaceObject (without a predefined space). It is possible to modify the default
 * space definition at any time. Space definition of pre-existing ASpaceObjects remains the same.
 * (no more shared pointer)
 */
class GSTLEARN_EXPORT ASpaceObject : public AStringable
{
public:
  ASpaceObject(const ASpaceSharedPtr& space = ASpaceSharedPtr());
  ASpaceObject(const ASpaceObject& r);
  ASpaceObject& operator= (const ASpaceObject& r);
  virtual ~ASpaceObject();

  /// AStringable interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

public:
  /// Accessor to the current object space context
  ASpaceSharedPtr getSpace() const { return _space; }
  /// Indicate if I am consistent with my current space context
  bool isConsistent() const { return isConsistent(_space); }

  void setSpace(ASpaceSharedPtr &&space) { _space = std::move(space); }

  /// Return unitary vector for the current space context
  VectorDouble getUnitaryVector() const;

  /// Indicate if I am consistent with the provided space
  bool isConsistent(const ASpaceSharedPtr& space) const;
  virtual bool isConsistent(const ASpace* space) const = 0;

  //////////////////////////////////////////////////////////
  /// Shortcuts to ASpace methods

  /// Return the number of dimension of the current space context
  unsigned int getNDim(int ispace = -1) const;
  /// Return the current space context origin coordinates
  const VectorDouble& getOrigin(int ispace = -1) const;

  /// Return the distance between two space points for the current space context
  double getDistance(const SpacePoint& p1,
                     const SpacePoint& p2,
                     int ispace = 0) const;

  /// Return all the distances (space composits) between two space points for the current space context
  VectorDouble getDistances(const SpacePoint& p1,
                            const SpacePoint& p2) const;

  /// Return the increment vector between two space points for the current space context
  VectorDouble getIncrement(const SpacePoint& p1,
                            const SpacePoint& p2,
                            int ispace = 0) const;
  void getIncrementInPlace(const SpacePoint& p1,
                           const SpacePoint& p2,
                           VectorDouble& ptemp,
                           int ispace = -1) const;

protected:
  /// Modify the Space dimension of an already created item (and create RN space)
  /// (To be used only during creation ... in particular when reading NF)
  void setNDim(int ndim);

protected:
  /// Current space context of the object
  ASpaceSharedPtr _space;

private:
  /// Dummy vector for the origin
  VectorDouble _dummy;
};

/// (Re)Defining the unique default global space
GSTLEARN_EXPORT void defineDefaultSpace(const ESpaceType& type,
                                        unsigned int ndim = 2,
                                        double param      = 0.);
/// Set the unique default global space from another one
GSTLEARN_EXPORT void setDefaultSpace(const ASpaceSharedPtr& space);

/// Return a clone of the unique default global space

GSTLEARN_EXPORT ESpaceType getDefaultSpaceType();
GSTLEARN_EXPORT int getDefaultSpaceDimension();
GSTLEARN_EXPORT const ASpace* getDefaultSpace();
GSTLEARN_EXPORT ASpaceSharedPtr getDefaultSpaceSh();

GSTLEARN_EXPORT bool isDefaultSpaceSphere();
