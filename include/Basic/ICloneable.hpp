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
#include <assert.h>

/**
 * Inherits from this interface to make your class cloneable.
 */
class GSTLEARN_EXPORT ICloneable
{
public:
  ICloneable() = default;
  virtual ~ICloneable() = default;

  virtual ICloneable* clone() const = 0;
};

// from https://stackoverflow.com/questions/65916601/clone-derived-class-from-base-class-pointer
template<typename Base, typename Derived>
class TCloneable: public Base
{
  Base* clone() const override
  {
    return new Derived(*static_cast<Derived*>(this));
  }
};
